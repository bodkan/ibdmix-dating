#!/bin/env bash

set -e

for chrom in $(seq 1 22); do

  echo "Processing chromosome $chrom..."

  vcf="/maps/projects/racimolab/data/mesoneo/20200629-impute.1000g.sg_freeze_v1/vcf/20200629-impute.1000g.sg_freeze_v1/${chrom}.neo.impute.1000g.vcf.gz"

  gt="gt_chr${chrom}.tsv"

  echo "Parsing VCF file '$vcf'..."

  # write a header line of the genotype table, including sample names
  printf "chrom\tpos\tref\talt\tpan_troglodytes\tAA\t" > $gt
  bcftools query -l $vcf | tr '\n' '\t' | sed 's/\t$/\n/' >> $gt

  # convert the full VCF file into a simple tab-separated table of GTs
  # after first subsetting to only biallelic SNPs
  bcftools view -m2 -M2 -v snps $vcf \
      | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/pan_troglodytes\t%INFO/AA[\t%GT]\n" \
      | awk 'BEGIN {FS=OFS="\t"}
             {
                 $5 = toupper($5); $6 = toupper($6);
		 split($6, aa_parts, "|"); $6=aa_parts[1];
		 if ($5 == $3) $5 = "0"; else if ($5 == $4) $5 = "1"; else $5 = "-";
	         if ($6 == $3) $6 = "0"; else if ($6 == $4) $6 = "1"; else $6 = "-";
		 print
	     }' \
      | sed 's#\.|\.#-#g' \
      | sed 's#\./\.#-#g' \
      >> $gt

  echo "Filtering for ancestry-informative markers..."

  Rscript ancestry_sites.R $gt

  rm $gt

done

zcat info_gt_chr{1..22}.tsv.gz \
  | awk 'NR == 1 && $2 == "pos" || NR > 1 && $2 != "pos"' \
  | bgzip > data/info_gt.tsv.gz
tabix -s1 -b2 -e2 data/info_gt.tsv.gz

rm info_*.tsv.gz
