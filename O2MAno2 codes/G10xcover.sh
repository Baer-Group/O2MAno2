#!/bin/sh

#SBATCH --job-name=G10xcover # Job name

#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail

#SBATCH --cpus-per-task=3 # Number of CPU cores per task

#SBATCH --ntasks=1 # Run a single task

#SBATCH --mem=10gb # Memory limit

#SBATCH --time=80:00:00 # Time limit hrs:min:sec

#SBATCH --output=G10xcover_%j.out # Standard output and error log

#SBATCH --account=baer --qos=baer



# Load the necessary module for samtools

module load samtools



# Constant reference genome length for C. elegans

genome_length=100286401



# Define the input directory where the .markduplicates.bam files are located

input_dir="/orange/baer/G_new"



# Output CSV file in the current working directory

output_file="G10x_coverage_results.csv"



# Write the header to the CSV file

echo "Sample,Total_Covered_Length,Breadth_of_Coverage(%),Total_Reads,Average_Depth" > $output_file



# Loop through all .markduplicates.bam files in the input directory

find "$input_dir" -type f -name "*.markduplicates.bam" | while read -r bam_file; do

  # Get the sample name by removing the .markduplicates.bam extension

  sample_name=$(basename "$bam_file" .markduplicates.bam)



  # Use samtools to calculate depth and store it in a temporary file

  samtools depth -a "$bam_file" > g_temp.txt



  # Calculate total covered length (number of bases with coverage >= 3)

  total_covered_length=$(awk '{ if ($3 >= 10) covered += 1 } END { print covered }' g_temp.txt)



  # Calculate total reads and average depth

  total_reads=$(awk '{sum += $3} END {print sum}' g_temp.txt)

  average_depth=$(awk -v total="$total_reads" -v genome="$genome_length" 'BEGIN { print (total / genome) }')



  # Calculate breadth of coverage as a fraction

  breadth_of_coverage=$(echo "scale=2; ($total_covered_length / $genome_length) * 100" | bc)



  # Append results to the CSV file

  echo "$sample_name,$total_covered_length,$breadth_of_coverage,$total_reads,$average_depth" >> $output_file



  # Clean up temporary file

  rm g_temp.txt



  echo "Processed: $sample_name"

done



echo "All coverage results saved to $output_file"


