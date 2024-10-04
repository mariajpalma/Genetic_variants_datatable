library(data.table)

### Chromosome 22
vcf_chr22 <- fread("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", 
                   skip = "#CHROM", 
                   na.strings = ".", 
                   select = c(1,2,4,5,sample(x = c(10:2513), size = 100, replace = FALSE)))

# Rename and reorganize columns
setnames(vcf_chr22, "#CHROM", "CHROM")
vcf_chr22[, ID := paste(CHROM, POS, sep = "_")]
setcolorder(vcf_chr22, c(colnames(vcf_chr22)[1:2], "ID", colnames(vcf_chr22)[3:(ncol(vcf_chr22) - 1)]))

### Step 1: Keep variant sites (remove fixed sites)
vcf_chr22 <- vcf_chr22[!is.na(ALT)]

### Step 2: Keep SNPs only (A, T, C, G in ALT)
vcf_chr22 <- vcf_chr22[grepl("^[ATCG]$", ALT)]
vcf_chr22 <- vcf_chr22[grepl("^[ATCG]$", REF)]

### Step 3: Filter variants with more than 5% missing data
genotype_cols <- names(vcf_chr22)[6:ncol(vcf_chr22)] # columns with genotype data

# Create a temporary column for counting missing genotypes
vcf_chr22[, missing_count := rowSums(.SD == ".", na.rm = TRUE), .SDcols = genotype_cols]

# Calculate total number of samples (excluding the fixed columns)
num_samples <- length(genotype_cols)

# Calculate the percentage of missing data
vcf_chr22[, missing_pct := missing_count / (num_samples * 2)]  # Each sample has two alleles

# Keep rows with less than 5% missing data
vcf_chr22 <- vcf_chr22[missing_pct < 0.05]

# Optionally remove the temporary columns
vcf_chr22[, `:=`(missing_count = NULL, missing_pct = NULL)]

### Step 4: Identify individuals with <= 5% missing data
missing_counts <- vcf_chr22[, lapply(.SD, function(x) sum(x == ".")), .SDcols = genotype_cols]

# Calculate total number of variants
total_variants <- nrow(vcf_chr22)

# Calculate percentage of missing data for each individual
missing_pct <- colSums(missing_counts) / total_variants

# Identify individuals to keep (<= 5% missing data)
individuals_to_keep <- names(missing_pct[missing_pct < 0.05])

# Filter the dataset to keep only relevant columns
info <- c("CHROM", "POS", "ID", "REF", "ALT")
info_ikeep <- c(info, individuals_to_keep)
vcf_chr22 <- vcf_chr22[, ..info_ikeep]

### Step 5: Count the number of alternate alleles per individual
# The genotype columns contain values like '0|0', '0|1', '1|1', etc.
# We want to count how many '1's (i.e., alternate alleles) appear for each individual.



# Split the genotype by "|" and count the number of '1's for each individual
vcf_chr22[, (genotype_cols) := lapply(.SD, function(gt) {
  # Split the genotype string (e.g., "0|1") and count the number of '1's
  sapply(strsplit(gt, "\\|"), function(alleles) sum(alleles == "1", na.rm = TRUE))
}), .SDcols = genotype_cols]

# View the first few rows with alternate allele counts
head(vcf_chr22)

### Read info about the samples
metadata <- fread("igsr_samples_1000G.tsv")
setnames(metadata, "Sample name", "Sample_name")
setnames(metadata, "Population code", "Population_code")
metadata <- metadata[, .(Sample_name, Population_code)]

# Melt the VCF data to long format
vcf_chr22_long <- melt(vcf_chr22, 
                 id.vars = c("CHROM", "POS", "ID", "REF", "ALT"),  # Columns to keep as identifiers
                 measure.vars = individuals_to_keep,                     # Genotype columns (samples)
                 variable.name = "Sample_name",                    # Name of the new 'sample' column
                 value.name = "Genotype")                          # Values will be in 'Genotype'

### Merge

vcf_chr22_long <- merge(vcf_chr22_long, metadata, by = "Sample_name", all.x = TRUE)

# Calculate the number of alternative alleles per variant and per population
variant_summary <- vcf_chr22_long[, .(
  Total_alt_alleles = sum(Genotype, na.rm = TRUE),
  Num_individuals = uniqueN(Sample_name)
  # Sum the Genotype values to get total alternative alleles for each variant
  ), by = .( ID, Population_code)]  # Group by variant and population

variant_summary[, alt_freq:= Total_alt_alleles/(Num_individuals*2)]


