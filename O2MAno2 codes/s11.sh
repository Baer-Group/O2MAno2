#!/bin/sh
#SBATCH --job-name=s11 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=9 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=30gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s11_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /blue/baer/m.rifat/G
### step 11: Apply 3x-coverage filtering:
module load gatk/
for i in *.vcf; do
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar VariantFiltration \
-R c_elegans.PRJNA13758.WS289.genomic.fa \
-O "${i/%.vcf/.3x.vcf}" \
--variant "$i" \
--genotype-filter-name "DP" \
--genotype-filter-expression "DP < 3.0"
done
for i in *.3x.vcf; do
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar SelectVariants \
-V "$i" \
--set-filtered-gt-to-nocall \
-O "${i/%.3x.vcf/.3x.nocall.vcf}"
done
for i in *.3x.nocall.vcf; do
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar SelectVariants \
-V "$i" \
-R c_elegans.PRJNA13758.WS289.genomic.fa \
-O "${i/%.3x.nocall.vcf/.3x.filtered.vcf}" \
--exclude-filtered true \
--exclude-non-variants true
done
