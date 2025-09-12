#!/usr/bin/env bash

echo -e "Sample1\tSample2\tCommonSites_3x" > common_3x_summary.txt



while IFS=$'\t' read -r bam1 bam2; do

  if [[ -z "$bam2" ]]; then

    echo -e "$(basename "$bam1")\tNA\t0" >> common_3x_summary.txt

    continue

  fi



  s1=${bam1%.bam}; s2=${bam2%.bam}

  shared=$(samtools depth -a "$bam1" "$bam2" | awk '$3>=3 && $4>=3' | wc -l)



  echo -e "${s1}\t${s2}\t${shared}" >> common_3x_summary.txt

done < pairs.txt


echo "Done – see common_3x_summary.txt"


