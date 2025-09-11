#! /bin/sh
#SBATCH --job-name=s5 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=9 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=80gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s5_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer 
cd /blue/baer/m.rifat/G
### step 5: Add read group information using picard:
module load picard/2.25.5/
for i in *.PP.bam; do
x=`basename $i .PP.bam`
picard AddOrReplaceReadGroups \
-I "$i" \
-O "${i/%.PP.bam/.PP.RG.bam}" \
-RGID BaerLines_WISeq \
-RGLB "$x" \
-RGPL illumina \
-RGPU BaerLines_WISeq \
-RGSM "$x"
done
