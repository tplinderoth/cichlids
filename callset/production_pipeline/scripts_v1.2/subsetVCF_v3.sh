#!/bin/bash

# subsetVCF_v3.sh [chr1 vcf file] [sample subset file] [reference fasta] [outfile prefix] [','-separated chromosome numbers to process]

if [ "$#" -lt 5 ]
then
	printf "\n%s\n\n" "subsetVCF_v3.sh [chr1 vcf file] [sample subset file] [reference fasta] [outfile prefix] [','-separated chromosomes numbers to process]"
	exit
fi

VCFTMP=$1
SUBFILE=$2
REF=$3
OUTPREFIX=$4
CHR=(${5//,/ })

MASKSTARS='/home/tl483/rds/rds-rd109-durbin-group/projects/cichlid/cichlid_g.vcf/main_set_2019-07-10/fAstCal1.2/GenotypeCorrected/depth_filter/scripts/maskStarGeno.pl'

for i in "${CHR[@]}"
do
	VCF=$( echo "$VCFTMP" | sed s/chr1/chr$i/ )
	OUTFILE="${OUTPREFIX}_chr${i}.vcf.gz"

	# subset samples and trims alternate alleles not in subset, mask genotypes with star alleles, trim star alleles,
	# left-align and normalize indels, recalculate INFO (allele frequencies, DP, excess heterozygosity, and number of individuals with data),
	# retain INFO and FORMAT data relevant only to subset, filter out nonvariants

	bcftools view -O v -a -S $SUBFILE -U -c 1:nref $VCF | $MASKSTARS | bcftools view -O u -a -U -c 1:nref | bcftools norm -O u -f $REF | bcftools +fill-tags -- -t 'AN,AC,AF,MAF,NS,ExcHet,DP=sum(DP)' \
| bcftools annotate -O u -x ^INFO/AN,INFO/AC,INFO/AF,INFO/MAF,INFO/NS,INFO/ExcHet,INFO/DP,^FORMAT/GT,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PL | bcftools view -O z -i 'MAF[*] > 0' > $OUTFILE &

done
wait

exit
