#!/bin/env bash

vcf="/maps/projects/racimolab/data/mesoneo/20200629-impute.1000g.sg_freeze_v1/vcf/20200629-impute.1000g.sg_freeze_v1/21.neo.impute.1000g.vcf.gz"

chrom="21"

gt="gt_chr${chrom}.tsv"

# write a header line of the genotype table, including sample names
printf "chrom\tpos\tref\talt\t" > $gt
bcftools query -l $vcf | tr '\n' '\t' | sed 's/\t$/\n/' >> $gt

# convert the full VCF file into a simple tab-separated table of GTs
# after first subsetting to only biallelic SNPs
bcftools view -m2 -M2 -v snps $vcf \
    | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n" \
    | sed 's#\.|\.#-#g' \
    | sed 's#\./\.#-#g' \
    >> $gt

Rscript ancestry_sites.R $gt

rm $gt
