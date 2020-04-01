#!/bin/bash

# extractDP.sh [bgzipped chr1 VCF] [output directory]

VCFTMP=$1
OUTDIR=$2
PLOTSCRIPT='/home/tl483/rds/rds-rd109-durbin-group/projects/cichlid/cichlid_g.vcf/main_set_2019-07-10/fAstCal1.2/GenotypeCorrected/depth_filter/scripts/plot_dp.R'
chr=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23)
header="CHROM\tPOS\tDP"
OUTDIR=$( echo "$OUTDIR" | sed 's/\/*$/\//' )

for i in "${chr[@]}"
do
	VCF=$( echo "$VCFTMP" | sed s/chr1/chr$i/g )
	OUT=$({ echo "$OUTDIR" ; echo "$VCF" | perl -ne 'print $1 if /([^\/]+)$/' | sed s/chr${i}.vcf.gz/chr${i}_DP.txt/ ; } | tr "\n" " " | sed 's/ //')
	OUT=$( echo "$OUT" | sed 's/_variants3_/_variants_/' )
	echo -e $header > $OUT

	bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' $VCF >> $OUT &

done
wait

echo "Extracting DP finished --> moving onto plotting"

$PLOTSCRIPT

exit
