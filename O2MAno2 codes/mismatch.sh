#!/usr/bin/env bash

set -euo pipefail

shopt -s nullglob



TYPES=(snp indel)

COVS=(3x 10x)



for TYPE in "${TYPES[@]}"; do

  for COV in "${COVS[@]}"; do



    ECSV="E.${TYPE}.${COV}_O1.csv"

    GCSV="G.${TYPE}.${COV}_O1.csv"

    OUT="mismatch.${TYPE}.${COV}.csv"



    echo ">>> Building cross‑sublines mismatches for $TYPE/$COV → $OUT"



    awk -F, -v OFS=, '

      ############################################################################

      # Phase 1: read the E‑CSV, capture REF/ALT & every E‑sample’s 1/1 state   #

      NR==FNR {

        # on header line, discover E‑columns (two naming patterns)

        if (FNR==1) {

          for(i=1; i<=NF; i++){

            h=$i; gsub(/^"|"$/, "", h)

            # match E###_ at start or _E###_ in middle

            if (match(h, /(^|_)E([0-9]{3})_/, m)) {

              Ecol[i]=h

              # grab last two digits of the three-digit ID

              Eid[i]=substr(m[2],2,2)

            }

          }

          next

        }

        # for each data row: store REF/ALT + E genotypes

        chr=$1; gsub(/^"|"$/, "", chr)

        pos=$2; gsub(/^"|"$/, "", pos)

        key=chr":"pos

        ref[key]=$3; gsub(/^"|"$/, "", ref[key])

        alt[key]=$4; gsub(/^"|"$/, "", alt[key])

        for(i in Ecol){

          ge=$i; gsub(/^"|"$/, "", ge)

          Eg[key,i]=ge

        }

        next

      }

      ############################################################################

      # Phase 2: read G‑CSV, compare each G’s 1/1 against all stored E’s         #

      FNR==1 {

        # header: discover G‑columns

        for(i=1; i<=NF; i++){

          h=$i; gsub(/^"|"$/, "", h)

          if (match(h, /(^|_)G([0-9]{3})_/, m)) {

            Gcol[i]=h

            Gid[i]=substr(m[2],2,2)

          }

        }

        # header of output

        print "Chromosome","Position","REF","ALT","E_sample","G_sample"

        next

      }

      {

        # for each data row in G‑CSV:

        chr=$1; gsub(/^"|"$/, "", chr)

        pos=$2; gsub(/^"|"$/, "", pos)

        key=chr":"pos

        if (!(key in ref)) next            # skip if E never saw it

        for(ei in Ecol){

          ge=Eg[key,ei]

          if (ge=="1/1"||ge=="1|1"){

            for(gi in Gcol){

              gg=$gi; gsub(/^"|"$/, "", gg)

              if (gg=="1/1"||gg=="1|1"){

                # only if last‑two‑digit IDs differ

                if (Eid[ei] != Gid[gi]) {

                  print chr, pos, ref[key], alt[key], Ecol[ei], Gcol[gi]

                }

              }

            }

          }

        }

      }

    ' "$ECSV" "$GCSV" > "$OUT"



    echo "    → done: $(wc -l <"$OUT") total mismatches"

    echo

  done

done


