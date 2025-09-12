#!/bin/sh
#SBATCH --job-name=s7 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=9 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=80gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s7_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /blue/baer/m.rifat/G
### step 7: Call variants per-sample using HaplotypeCaller (in BP_RESOLUTION mode):
module load gatk/
module load samtools/
module load picard/
samtools faidx c_elegans.PRJNA13758.WS289.genomic.fa
picard CreateSequenceDictionary -R c_elegans.PRJNA13758.WS289.genomic.fa -O c_elegans.PRJNA13758.WS289.genomic.dict
for i in *.markduplicates.bam; do
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller \
-R c_elegans.PRJNA13758.WS289.genomic.fa \
-I "$i" \
-O "${i/%.markduplicates.bam/.bp.g.vcf.gz}" \
-ERC BP_RESOLUTION
done
