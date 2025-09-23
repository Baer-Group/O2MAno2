Directory Description:
Chromosome_distribution: Contains the location (Chromosome+Position) of all the mutations.
O2MAno2 Codes: Contains all the codes used for genome assembly, variant calling, annotation and mutation rate & spectra. .
Pymc_scripts: Contains all the codes used to analyze fitness data.
dn_ds: Contains dN/dS analysis dataset and code. 
merge_EG: Contains total counts of mutation per sample at different coverage (3X and 10x) (SNP/Indel) and regions (whole genome, mononucleotide repeat and non repeat)
Workflow:

1.	s2.sh Gunzip to decompress fq.gz files.
 fastp to trim the adapters and the low-quality reads from the sample genome files.
2.	samtools faidx to create. fa index of the reference genome WS292 (command line)
3.	picard createsequencedictionary for creating dictionary for the reference genome (.dict file) (command line)
4.	s3.sh bowtie2 for aligning the reads with the reference genome.
samtools view to convert sam(Sequence Alignment Map) files into bam(Binary Alignment Map).
samtools sort to sort bam files in the order of the chromosome and positions.
bamtools filter -isProperPair true to make sure all the read pairs are aligned and in a correct order. 
5.	s5.sh picard AddOrReplaceReadGroups  add group information to the files. 
6.	s6.sh MarkDuplicates command in gatk to remove duplicate reads.
7.	sq.sh fastqc and multiqc to check the read quality.
8.	s7.sh HaplotypeCaller in gatk to call the variants.
9.	s8.sh  GenomicsDBImport in gatk to consolidate all the vcf files in a database.
10.	s9.sh GenotypeGVCFs in gatk to call variants jointly and generate vcf files.
11.	s10.sh SelectVariants in gatk to extract SNPs and INDELs.
12.	s11.sh VariantFiltration in gatk to apply 3X and 10X coverage filtering.
      SelectVariants to fileter no call variants. 
13.	s12.sh bcftools view max allele to convert into bi allelic vcf files 
14.	mutation.sh bcftools view count to filter out the rows other than the mutation rows.
15.	pair.sh List all the subline pairs in a text file named “pair.txt”. 
16.	combine2.sh and cov.sh  Calculate the number of base pairs covered in both sublines using the text file “pair.txt”.. Record the coverage data in two text files named “common_3x_summary.txt” and “common_10x_summary.txt”
17.	vcf_to_csv.R to convert all the files into csv files. (package used: dplyr, vcfR)
18.	tandem repeat finder (trf) 2 5 5 75 10 0 2 -d -m -h to identify repeat regions of the reference genome WS292. (command line)
19.	mono75.R, to extract the mono and dinucleotide repeat regions in the reference region and convert to .csv files. (Total mono repeats, 6786116) (package used: dplyr) (use dat file from previous step)
20.	mono.R split all the files into two groups: mono and non mono (package used: dplyr, vcfR). 
21.	Order.R, to split all the files into two groups again, o1 and o2) (package used: dplyr)
22.	indel.R and snp.R, to distribute the mutations into indel and snp spectra. Calculated total for each sample (package used: dplyr). 

23.	3Xcover.sh, 10Xcover.sh, G3xcover.sh, G10xcover.sh using these four-custom made bash script to calculate the percentage of the total length covered by the MA genomes. 
24.	cat.R to add coverage information and delete ancestors. 
25.	merge_new.R Combines the snp and indel files. 
26.	merge_EG.R Combines the E and G sublines in a single file. 
27.	Indel_spectra.R Plots insertion and deletion spectra seperately. 
28.	Irate_spectra.R Plots SNV and Indel spectra seperately. 
29.	
30.	Dummy_reference.R . Randomly inserted SNV and indels to generate a pseudo reference genome. (package used: GenomicRanges,IRanges,Biostrings)
31.	GATK: FastaAlternateReferenceMaker. Pseudo-Reference genome created. (command line). 
32.	fn.R To calculate False negatives and failure to recall. 
33.	sim_plain.R To calculate the point estimates of the rates per generation/days with and without FP correction. 
34.	sim_from_o1.R Testing the hypothesis of uniform mean.
sim_from_o1_v2.R Testing the hypothesis of uniform mean(Updated version).
35.	sd_seperate_o1o2.R Testing the hypothesis of uniform variance in o1 and o2. 
36.	mismatch.sh . Diagnose the mismatched o1 lines. 
37.	snpeff . Annotate the mutations, and their putative impacts. (Download and building database command changed)
38.	annotation.R and add_line.R . Process the annotated vcf files and convert to csv file.
39.	ann.sh  . Generate tsv files from annotated vcf files.  
40.	annotation_v2.R . Annotate genes from tsv files and update the annotation CSV file.
41.	mutation_distribution_v2.R . Distribution of the mutations across the chromosomes(package used: readr). 
42.	dn_ds.R . Used to convert o1 o2 csv to vcf to annotate and calculate dn/ds(package used: stringr, forcats). 
43.	tsv.sh . To generate annotated tsv files from the previous steps.  
44.	working_command.txt . Use the command from this text directly and change manually for two files in case tsv.sh doesn’t work.  
45.	dn_ds_v2.R . To generate the dN/dS table from annotated tsv files. 



