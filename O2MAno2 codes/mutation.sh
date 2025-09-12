#!/bin/bash

module load bcftools

# Loop through all .vcf files in the current directory

for vcf_file in *filtered.biallelic.vcf

do

    # Define the output file name by replacing the original filename's extension

    output_file1="${vcf_file%filtered.biallelic.vcf}mutation.vcf"
    output_file2="${vcf_file%filtered.biallelic.vcf}reversion.vcf"

    

    # Run bcftools command on each file

    bcftools view -i 'COUNT(GT="alt")=1 && COUNT(GT="het")=0 && COUNT(GT="miss")<=10' "$vcf_file" -Oz -o "$output_file1"
    bcftools view -i 'COUNT(GT="ref")=1 && COUNT(GT="het")=0 && COUNT(GT="miss")<=10' "$vcf_file" -Oz -o "$output_file2"


    echo "Processed $vcf_file -> $output_file1"
    echo "Processed $vcf_file -> $output_file2"

done


