#!/usr/bin/env bash
module load bcftools
# generate pairs.txt of E4* vs G2* markduplicates.bam

out="pairs.txt"

: > "$out"   # truncate or create



# loop over every E4…markduplicates.bam

for e in *E4*.markduplicates.bam; do

  # extract the E-ID (e.g. "E400", "E412", etc.)

  eid=$(basename "$e" | grep -oE 'E4[0-9]+' )

  # map to the corresponding G-ID (E4NNN → G2NNN)

  gid=${eid/E4/G2}



  # try to find a G-file containing that gid

  g=$(ls *"${gid}"*.markduplicates.bam 2>/dev/null | head -n1)



  if [[ -n "$g" ]]; then

    echo -e "${e}\t${g}" >> "$out"

  else

    # leave second column blank (or you could write "NONE")

    echo -e "${e}\t" >> "$out"

  fi

done



echo "Done — see $out for your pairs."


