#!/usr/bin/bash
## Calculate summary statistics
VCF=$1
OUTDIR=$2
RTG=~/sw/rtg-core-non-commercial-3.9.1/rtg

## Index VCF
tabix -f ${VCF}

## Get whole genome variant stats
# mkdir ${OUTDIR} ## dones in R only for testing script
${RTG} vcfstats ${VCF} > ${OUTDIR}/chr1-22.stats.txt


## Get chromsome level stats
for i in {1..22}; do
    tabix -h ${VCF} ${i} > ${OUTDIR}/chr${i}.vcf
    ${RTG} vcfstats ${OUTDIR}/chr${i}.vcf > ${OUTDIR}/chr${i}_stats.txt
    rm ${OUTDIR}/chr${i}.vcf
done