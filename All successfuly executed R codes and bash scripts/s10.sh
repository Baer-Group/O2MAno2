#!/bin/sh
#SBATCH --job-name=s10 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=9 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=30gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s10_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /blue/baer/m.rifat/G
### step 10: Extract SNPs & INDELs:
module load gatk/
for i in *.vcf; do
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar SelectVariants \
-R c_elegans.PRJNA13758.WS289.genomic.fa \
-V "$i" \
 --select-type-to-include SNP \
 -O "${i/%.vcf/.snp.vcf}"
done
for i in *.vcf; do
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar SelectVariants \
-R c_elegans.PRJNA13758.WS289.genomic.fa \
-V "$i" \
--select-type-to-include INDEL \
-O "${i/%.vcf/.indel.vcf}"
done
