# EDIT: point to your reference

REF=/blue/baer/m.rifat/STR/c_elegans.PRJNA13758.WS292.genomic.fa



for VCF in E.bp.GVCFs.snp.3x.mutation.ann.vcf G.bp.GVCFs.snp.3x.mutation.ann.vcf E.bp.GVCFs.indel.3x.mutation.ann.vcf G.bp.GVCFs.indel.3x.mutation.ann.vcf

do

  base=${VCF%.vcf}



  # 1) Normalize & split multiallelics to one ALT per row

  bcftools norm -m -both -f "$REF" -Oz -o ${base}.split.vcf.gz "$VCF"

  tabix -p vcf ${base}.split.vcf.gz



  # 2) Emit CHROM POS REF ALT and the raw ANN string;

  #    then split ANN (comma) so each ANN record becomes a separate line

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' ${base}.split.vcf.gz | awk -F'\t' 'BEGIN{OFS="\t"}{

      n=split($5,a,",");

      for(i=1;i<=n;i++) print $1,$2,$3,$4,a[i]

    }' > ${base}.ANN_expanded.tsv

done


