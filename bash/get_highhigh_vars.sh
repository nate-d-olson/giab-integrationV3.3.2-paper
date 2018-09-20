#!/usr/bin/bash
BED=${1}
VCF=${2}
OUTDIR=${3}

## Get the intersection of high confidence vcf and bed files
bedtools intersect -header -a ${VCF} -b ${BED} > ${OUTDIR}/highhigh.vcf
bgzip ${OUTDIR}/highhigh.vcf
tabix -f ${OUTDIR}/highhigh.vcf.gz