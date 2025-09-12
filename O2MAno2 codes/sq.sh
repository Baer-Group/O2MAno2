#!/bin/sh
#SBATCH --job-name=sq # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=2 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=sq_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /orange/baer/CB_MA_seq_B/E
module load fastqc/
fastqc *markduplicates.bam
