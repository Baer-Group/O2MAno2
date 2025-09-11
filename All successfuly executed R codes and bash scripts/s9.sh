#!/bin/sh
#SBATCH --job-name=s9 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=9 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=80gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s9_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /blue/baer/m.rifat/G
### step 9: Call variants jointly using GenotypeGVCFs in GATK:
module load gatk/
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar GenotypeGVCFs -R c_elegans.PRJNA13758.WS289.genomic.fa -V gendb://G -O G.bp.GVCFs.vcf
