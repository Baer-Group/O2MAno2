library(vcfR)

library(dplyr)



# Define a list of VCF file paths

vcf_files <- list.files(pattern = "*.vcf")  # Adjust the pattern if your files have a different extension



# Loop over each VCF file

for (vcf_file in vcf_files) {

  

  # Read each VCF file

  vcf <- read.vcfR(vcf_file)

  

  # Extract the chromosome, position, reference, and alternate allele information

  chrom <- vcf@fix[, "CHROM"]

  pos <- vcf@fix[, "POS"]

  ref <- vcf@fix[, "REF"]

  alt <- vcf@fix[, "ALT"]

  

  # Extract genotype information

  gt <- extract.gt(vcf, element = "GT")

  

  # Create a data frame with the extracted information

  df <- data.frame(

    Chromosome = chrom,

    Position = pos,

    REF = ref,

    ALT = alt,

    gt,

    stringsAsFactors = FALSE  # Ensure genotypes are treated as characters

  )

  

  # Generate a corresponding CSV file name for each VCF file

  output_csv <- sub(".vcf$", ".rows.csv", vcf_file)

  

  # Save the summary table to a CSV file

  write.csv(df, output_csv, row.names = FALSE)

  

  # Print a message indicating the file has been processed

  cat("Processed:", vcf_file, "->", output_csv, "\n")

}


