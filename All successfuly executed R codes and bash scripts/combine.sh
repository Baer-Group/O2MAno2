#!/usr/bin/env bash

shopt -s nullglob



TYPES=(snp indel)

COVS=(3x 10x)



for TYPE in "${TYPES[@]}"; do

  for COV in "${COVS[@]}"; do



    # 1) all E→G missing beds for this TYPE/COV

    cat *."${TYPE}"."${COV}".missing_in_G*.bed | sort -k1,1 -k2,2n | uniq \

      > combined."${TYPE}"."${COV}".E_missing.bed



    # 2) all G→E missing beds for this TYPE/COV

    cat *."${TYPE}"."${COV}".missing_in_E*.bed | sort -k1,1 -k2,2n | uniq \

      > combined."${TYPE}"."${COV}".G_missing.bed



    echo "Written:"

    echo "  combined.${TYPE}.${COV}.E_missing.bed"

    echo "  combined.${TYPE}.${COV}.G_missing.bed"

    echo



  done

done


