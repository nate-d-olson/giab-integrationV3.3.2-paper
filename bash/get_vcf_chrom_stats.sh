#!/usr/bin/bash
## Calculate summary statistics
VCF=${1}
OUTDIR=${2}

## Get whole genome variant stats
   rtg vcfstats ${VCF} > ${OUTDIR}/chr1-22.stats.txt
   
## Get chromsome level stats
for i in {1::22}; do
    tabix -h ${VCF} ${i} > ${OUTDIR}chr${i}.vcf
    rtg vcfstats chr${j}.vcf.gz > ${OUTDIR}/chr${i}_stats.txt
    rm ${OUTDIR}/chr${j}.vcf
done