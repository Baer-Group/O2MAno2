#!/usr/bin/env bash

set -euo pipefail



# ====== CONFIG (edit these paths if needed) ======

REF="/blue/baer/m.rifat/STR/c_elegans.PRJNA13758.WS292.genomic.fa"

VCF_O1="O1.snps.ann.vcf.gz"

VCF_O2="O2.snps.ann.vcf.gz"

OUT="ann_tsv_out"



# Optional: two-column contig rename map (leave empty if not needed)

#   chrI  I

#   chrII II

CHR_MAP="${CHR_MAP:-}"  # e.g., export CHR_MAP=chr.map



# ====== deps ======

module load samtools >/dev/null 2>&1 || true

module load bcftools >/dev/null 2>&1 || true



mkdir -p "$OUT"



echo "[*] Indexing reference"

samtools faidx "$REF"



normalize_vcf () {

  local in="$1" out="$2"

  echo "[*] Normalizing $in -> $out"

  if [[ -n "${CHR_MAP}" && -s "${CHR_MAP}" ]]; then

    bcftools annotate --rename-chrs "${CHR_MAP}" -Ou "$in" | bcftools reheader -f "${REF}.fai" -Ou | bcftools view -Oz -o "$out"

  else

    bcftools reheader -f "${REF}.fai" "$in" | bcftools view -Oz -o "$out"

  fi

  tabix -f -p vcf "$out"

}



# 1) Reheader → BGZF → index

normalize_vcf "$VCF_O1" "${OUT}/O1.snps.ann.hdr.vcf.gz"

normalize_vcf "$VCF_O2" "${OUT}/O2.snps.ann.hdr.vcf.gz"



# 2) Extract ANN to TSV

echo "[*] Extracting ANN to TSV"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "${OUT}/O1.snps.ann.hdr.vcf.gz" > "${OUT}/O1.ANN.tsv"



bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' "${OUT}/O2.snps.ann.hdr.vcf.gz" > "${OUT}/O2.ANN.tsv"



echo "[✓] Done. TSVs at: ${OUT}/O1.ANN.tsv  and  ${OUT}/O2.ANN.tsv"


