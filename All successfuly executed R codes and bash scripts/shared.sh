#!/usr/bin/env bash

set -euo pipefail

shopt -s nullglob



PAIRS="pairs.txt"

TYPES=(snp indel)

COVS=(3x 10x)

E_PREFIX="E.bp.GVCFs"

G_PREFIX="G.bp.GVCFs"



while read -r E_BAM G_BAM; do

  # Extract the 4‑digit IDs (e.g. E400, G200)

  E_ID=$(basename "$E_BAM" | grep -oE '[EG][0-9]{3}')

  G_ID=$(basename "$G_BAM" | grep -oE '[EG][0-9]{3}')



  for TYPE in "${TYPES[@]}"; do

    for COV in "${COVS[@]}"; do

      thr=${COV%x}   # drop the 'x' → 3 or 10



      E_VCF="${E_PREFIX}.${TYPE}.${COV}.mutation.vcf.gz"

      G_VCF="${G_PREFIX}.${TYPE}.${COV}.mutation.vcf.gz"



      # Pull out the exact sample name from the VCF header

      SAMPLE_E=$(bcftools query -l "$E_VCF" | grep -E '(^|_)'"${E_ID}"'_')

      SAMPLE_G=$(bcftools query -l "$G_VCF" | grep -E '(^|_)'"${G_ID}"'_')



      # E → missing in G

      OUT1="${E_ID}.${TYPE}.${COV}.missing_in_${G_ID}.bed"; bcftools query -s "$SAMPLE_E" -i 'GT="1/1"||GT="1|1"' -f '%CHROM\t%POS0\t%POS\n' "$E_VCF" | samtools depth -aa -b - "$G_BAM" | awk -v t=$thr '$3 < t {print $1"\t"$2-1"\t"$2}' > "$OUT1"




      # G → missing in E

      OUT2="${G_ID}.${TYPE}.${COV}.missing_in_${E_ID}.bed"; bcftools query -s "$SAMPLE_G" -i 'GT="1/1"||GT="1|1"' -f '%CHROM\t%POS0\t%POS\n' "$G_VCF" | samtools depth -aa -b - "$E_BAM" | awk -v t=$thr '$3 < t {print $1"\t"$2-1"\t"$2}' > "$OUT2"





      # Report counts

      echo "Pair $E_ID→$G_ID ($TYPE,$COV):"

      echo "  $OUT1 → $(wc -l < "$OUT1") sites"

      echo "  $OUT2 → $(wc -l < "$OUT2") sites"

      echo

    done

  done

done < "$PAIRS"


