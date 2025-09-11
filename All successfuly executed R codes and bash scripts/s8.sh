#!/bin/sh
#SBATCH --job-name=s8 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=9 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=80gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=s8_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /blue/baer/m.rifat/G
### step 8: consolidating the vcf.gz files using GenomicsDBImport:
module load gatk/
java -Xmx8g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar GenomicsDBImport \
-V G-1_POOLRET72_S213_L002.bp.g.vcf.gz \
-V G-2_POOLRET73_S469_L003.bp.g.vcf.gz \
-V G-4_POOLRET73_S411_L003.bp.g.vcf.gz \
-V G-5_POOLRET73_S412_L003.bp.g.vcf.gz \
-V G201_POOLRET73_S509_L003.bp.g.vcf.gz \
-V G202_POOLRET73_S446_L003.bp.g.vcf.gz \
-V G203_POOLRET72_S273_L002.bp.g.vcf.gz \
-V G204_POOLRET73_S504_L003.bp.g.vcf.gz \
-V G206_POOLRET73_S465_L003.bp.g.vcf.gz \
-V G209_POOLRET72_S352_L002.bp.g.vcf.gz \
-V G210_POOLRET73_S494_L003.bp.g.vcf.gz \
-V G211_POOLRET73_S395_L003.bp.g.vcf.gz \
-V G212_POOLRET73_S396_L003.bp.g.vcf.gz \
-V G213_POOLRET72_S193_L002.bp.g.vcf.gz \
-V G214_POOLRET72_S332_L002.bp.g.vcf.gz \
-V G215_POOLRET72_S264_L002.bp.g.vcf.gz \
-V G216_POOLRET73_S502_L003.bp.g.vcf.gz \
-V G217_POOLRET73_S397_L003.bp.g.vcf.gz \
-V G218_POOLRET73_S521_L003.bp.g.vcf.gz \
-V G219_POOLRET73_S574_L003.bp.g.vcf.gz \
-V G221_POOLRET72_S308_L002.bp.g.vcf.gz \
-V G222_POOLRET72_S297_L002.bp.g.vcf.gz \
-V G223_POOLRET73_S451_L003.bp.g.vcf.gz \
-V G224_POOLRET73_S441_L003.bp.g.vcf.gz \
-V G226_POOLRET72_S258_L002.bp.g.vcf.gz \
-V G227_POOLRET72_S377_L002.bp.g.vcf.gz \
-V G229_POOLRET73_S539_L003.bp.g.vcf.gz \
-V G230_POOLRET73_S458_L003.bp.g.vcf.gz \
-V G231_POOLRET73_S457_L003.bp.g.vcf.gz \
-V G232_POOLRET73_S558_L003.bp.g.vcf.gz \
-V G233_POOLRET72_S242_L002.bp.g.vcf.gz \
-V G234_POOLRET72_S215_L002.bp.g.vcf.gz \
-V G235_POOLRET73_S570_L003.bp.g.vcf.gz \
-V G237_POOLRET72_S307_L002.bp.g.vcf.gz \
-V G239_POOLRET73_S444_L003.bp.g.vcf.gz \
-V G241_POOLRET73_S452_L003.bp.g.vcf.gz \
-V G242_POOLRET72_S286_L002.bp.g.vcf.gz \
-V G244_POOLRET73_S408_L003.bp.g.vcf.gz \
-V G245_POOLRET73_S519_L003.bp.g.vcf.gz \
-V G247_POOLRET73_S530_L003.bp.g.vcf.gz \
-V G248_POOLRET73_S387_L003.bp.g.vcf.gz \
-V G249_POOLRET72_S198_L002.bp.g.vcf.gz \
-V G252_POOLRET73_S481_L003.bp.g.vcf.gz \
-V G253_POOLRET73_S431_L003.bp.g.vcf.gz \
-V G254_POOLRET73_S447_L003.bp.g.vcf.gz \
-V G255_POOLRET72_S237_L002.bp.g.vcf.gz \
-V G257_POOLRET73_S540_L003.bp.g.vcf.gz \
-V G258_POOLRET73_S473_L003.bp.g.vcf.gz \
-V G259_POOLRET73_S428_L003.bp.g.vcf.gz \
-V G261_POOLRET72_S254_L002.bp.g.vcf.gz \
-V G262_POOLRET73_S394_L003.bp.g.vcf.gz \
-V G264_POOLRET73_S470_L003.bp.g.vcf.gz \
-V G265_POOLRET72_S298_L002.bp.g.vcf.gz \
-V G266_POOLRET73_S407_L003.bp.g.vcf.gz \
-V G268_POOLRET72_S342_L002.bp.g.vcf.gz \
-V G269_POOLRET73_S573_L003.bp.g.vcf.gz \
-V G270_POOLRET72_S345_L002.bp.g.vcf.gz \
-V G271_POOLRET73_S413_L003.bp.g.vcf.gz \
-V G273_POOLRET73_S525_L003.bp.g.vcf.gz \
-V G274_POOLRET72_S271_L002.bp.g.vcf.gz \
-V G277_POOLRET72_S364_L002.bp.g.vcf.gz \
-V G279_POOLRET72_S305_L002.bp.g.vcf.gz \
-V G281_POOLRET73_S422_L003.bp.g.vcf.gz \
-V G282_POOLRET72_S248_L002.bp.g.vcf.gz \
-V G283_POOLRET73_S427_L003.bp.g.vcf.gz \
-V G284_POOLRET73_S463_L003.bp.g.vcf.gz \
-V G285_POOLRET72_S325_L002.bp.g.vcf.gz \
-V G286_POOLRET73_S548_L003.bp.g.vcf.gz \
-V G287_POOLRET73_S549_L003.bp.g.vcf.gz \
-V G289_POOLRET73_S537_L003.bp.g.vcf.gz \
-V G290_POOLRET73_S426_L003.bp.g.vcf.gz \
-V G295_POOLRET73_S434_L003.bp.g.vcf.gz \
-V G296_POOLRET73_S386_L003.bp.g.vcf.gz \
-V G297_POOLRET72_S331_L002.bp.g.vcf.gz \
-V G298_POOLRET72_S306_L002.bp.g.vcf.gz \
-V G299_POOLRET72_S349_L002.bp.g.vcf.gz \
-V H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_a.bp.g.vcf.gz \
-V H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_b.bp.g.vcf.gz \
-V H7GWFDSX7_n01_EG8072_BaerAncestor_POOLRET74_c.bp.g.vcf.gz \
--genomicsdb-workspace-path G \
--intervals I --intervals II --intervals III --intervals IV --intervals V --intervals X --intervals MtDNA
