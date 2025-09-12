#Set current folder of the R script as the working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()
library(vcfR)
library(dplyr)
library(tidyr)
##bash code
### download database
## snpeff -v  -download WBcel235.105 -dataDir /blue/baer/m.rifat/MA_lines/snpeff_data
### snpeff -dataDir /blue/baer/m.rifat/MA_lines/snpeff_data WBcel235.105 E.bp.GVCFs.snp.3x.mutation.vcf > E.bp.GVCFs.snp.3x.mutation.ann.vcf
# Load the annotated VCF file
vcf <- read.vcfR("G.bp.GVCFs.indel.3x.mutation.ann.vcf")

# Extract the fixed data, which includes CHROM, POS, REF, ALT, etc.
vcf_data <- as.data.frame(vcf@fix)

# Extract the INFO field from the VCF data
vcf_data$INFO <- vcf@fix[,"INFO"]

# Filter the rows where INFO contains the ANN field
vcf_data_with_ann <- vcf_data %>%
  filter(grepl("ANN=", INFO))

# Separate the ANN field into individual annotations
vcf_data_with_ann <- vcf_data_with_ann %>%
  separate(INFO, into = c("Allele", "Annotation", "Putative_impact", "Gene_name", "Gene_id", 
                          "Feature_type", "Feature_id", "Transcript_biotype"), 
           sep = "\\|", extra = "drop", fill = "right")

# Select relevant columns, such as chromosome, position, gene name, etc.
gene_list <- vcf_data_with_ann %>%
  select(CHROM, POS, REF, ALT, Gene_name, Putative_impact, Annotation)

# View the first few rows to verify
print(head(gene_list))

# Write the final gene list to a CSV file
write.csv(gene_list, "genes_from_indel_G.csv", row.names = FALSE)
