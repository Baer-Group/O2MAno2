#!/usr/bin/env bash

set -euo pipefail

shopt -s nullglob



TYPES=(snp indel)

COVS=(3x 10x)



for TYPE in "${TYPES[@]}"; do

  for COV in "${COVS[@]}"; do



    #–– combine all E→G beds for this TYPE/COV

    beds=( *."${TYPE}"."${COV}".missing_in_G*.bed )

    if (( ${#beds[@]} == 0 )); then

      echo "⚠️  No E→G BEDs for ${TYPE}/${COV}"

    else

      cat "${beds[@]}" | sort -k1,1 -k2,2n | uniq > combined."${TYPE}"."${COV}".E_missing.bed

      echo "✔️  combined.${TYPE}.${COV}.E_missing.bed → $(wc -l < combined.${TYPE}.${COV}.E_missing.bed) lines"

    fi



    #–– combine all G→E beds for this TYPE/COV

    beds=( *."${TYPE}"."${COV}".missing_in_E*.bed )

    if (( ${#beds[@]} == 0 )); then

      echo "⚠️  No G→E BEDs for ${TYPE}/${COV}"

    else

      cat "${beds[@]}" | sort -k1,1 -k2,2n | uniq > combined."${TYPE}"."${COV}".G_missing.bed

      echo "✔️  combined.${TYPE}.${COV}.G_missing.bed → $(wc -l < combined.${TYPE}.${COV}.G_missing.bed) lines"

    fi



    echo

  done

done


