#!/bin/sh
#SBATCH --job-name=sa # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=80gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=sa_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /blue/baer/m.rifat/E
module load bowtie2/
bowtie2 -x ws289 -U H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_a_trimmed.fastq --phred33 -p 8 --very-sensitive-local -S H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_a.sam
bowtie2 -x ws289 -U H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_b_trimmed.fastq --phred33 -p 8 --very-sensitive-local -S H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_b.sam
bowtie2 -x ws289 -U H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_c_trimmed.fastq --phred33 -p 8 --very-sensitive-local -S H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_c.sam
