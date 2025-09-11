#!/bin/sh

#SBATCH --job-name=vcf_to_csv # Job name

#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail

#SBATCH --cpus-per-task=2 # Number of CPU cores per task

#SBATCH --ntasks=1 # Run a single task

#SBATCH --mem=30gb # Memory limit

#SBATCH --time=120:00:00 # Time limit hrs:min:sec

#SBATCH --output=vcf_to_csv_%j.out # Standard output and error log

#SBATCH --account=baer --qos=baer



cd /blue/baer/m.rifat/MA_lines/mono



#### Modules ####

module load R


# Run the R script to process all VCF files and convert them to CSV files

Rscript vcf_to_csv.R


