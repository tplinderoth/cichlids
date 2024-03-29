Malawi Cichlid Callset v1.3
========================

* **2198 individuals** 
* **47 genera** 
* **255 species (unique taxa)** 
___________________________

All-sites VCF location:<br/>
/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/allsites_vcf3/annotated

All-sites, sites-only VCF (same as all-sites VCF but without genotype information):<br/>
/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/allsites_vcf3/annotated/sites_only

All-variants VCF location:<br/>
/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/all_variants3

Biallelic VCF location:<br/>
/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3

Phased SNP VCF (all biallelic SNPs):<br/>
/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/phase/all_biallelic

Phased INFO/PASS SNP VCF location:<br/>
/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/biallelic3/phase/pass_only

### Malawi cichlids v1.3 summary

* all-sites VCF: Number of site entries (lines discounting the header) in all-sites VCF. Note that the same position may occur on multiple, adjacent lines.
* variants: Number of SNPs + indels
* SNPs: Number of SNPs
* insertions: Number of insertions
* deletions: Number of deletions
* SNPs: Number of single nucleotide changes
* SNP_sites: Number of sites with at least one single nucleotide change
* biallelic_SNPs: Number of biallelic SNP sites
* deletion_SNPs: Number of SNPs within deletions
* bi_deletion_SNPs: Number of biallelic SNPs within deletions
* biallelic_VCF: Number of sites in biallelic VCF (biallelic SNPs outside of deletions and with INFO/AF in interval (0,1))
* ancestral_annotations: Number of lines with an inferred ancestral allele
* QC_pass: Number of sites passing quality controls (see VCF headers, [annotation_summary](./annotation_summary.txt), and [chrM_annotation_summary](./chrM_annotation_summary.txt) for a detailed breakdown)

| chromosome | all-sites_VCF | variants  | insertions | deletions | SNPs      | SNP_sites | biallelic_SNPs | deletion_SNPs | bi_deletion_SNPs | biallelic_VCF | ancestral_annotations | QC_pass           |
|:-----------|:--------------|:----------|:-----------|:----------|:----------|:----------|:---------------|:--------------|:-----------------|:--------------|:----------------------|:------------------|
chr1         | 42665030      | 9657841   | 814475     | 1087556   | 7755810   | 7199772   | 6661102        | 1276517       | 1152008          | 5508718       | 35085948 (82.2%)      | 34447433 (80.7%)  |
chr2         | 39526156      | 8951247   | 698358     | 980491    | 7272398   | 6741746   | 6227568        | 1181295       | 1064937          | 5162263       | 31383743 (79.4%)      | 29930724 (75.7%)  |
chr3         | 53504991      | 11932209  | 865768     | 1204022   | 9862419   | 9086864   | 8337698        | 1474662       | 1320029          | 7015781       | 34485129 (64.5%)      | 27056333 (50.6%)  |
chr4         | 34306062      | 7891228   | 641469     | 876681    | 6373078   | 5896862   | 5435956        | 1040949       | 935693           | 4499767       | 27273014 (79.5%)      | 25101311 (73.2%)  |
chr5         | 40137524      | 9374000   | 783899     | 1054328   | 7535773   | 6981064   | 6444073        | 1267460       | 1141049          | 5302706       | 33063765 (82.4%)      | 32076363 (79.9%)  |
chr6         | 42902476      | 9872095   | 835714     | 1113675   | 7922706   | 7340156   | 6776061        | 1320879       | 1190071          | 5585504       | 34641840 (80.7%)      | 32792210 (76.4%)  |
chr7         | 69859057      | 15843551  | 1337685    | 1796985   | 12708881  | 11788621  | 10897972       | 2095724       | 1888453          | 9008537       | 56727768 (81.2%)      | 54641431 (78.2%)  |
chr8         | 26818740      | 6554091   | 526000     | 734337    | 5293754   | 4884841   | 4489320        | 928647        | 832972           | 3656197       | 22459177 (83.7%)      | 21628032 (80.6%)  |
chr9         | 36968958      | 8763579   | 683499     | 936314    | 7143766   | 6591082   | 6056560        | 1193244       | 1069788          | 4986130       | 27647342 (74.8%)      | 25548766 (69.1%)  |
chr10        | 35164388      | 8085216   | 667169     | 909650    | 6508397   | 6027053   | 5561137        | 1099998       | 989808           | 4570908       | 28620600 (81.4%)      | 27233436 (77.4%)  |
chr11        | 37872216      | 8521725   | 727663     | 954448    | 6839614   | 6350853   | 5877439        | 1102096       | 995037           | 4882037       | 31603157 (83.4%)      | 30089945 (79.5%)  |
chr12        | 40056780      | 9175348   | 782181     | 1031723   | 7361444   | 6815639   | 6287159        | 1213492       | 1090474          | 5196285       | 33109659 (82.7%)      | 31331385 (78.2%)  |
chr13        | 34331999      | 7605569   | 640176     | 858342    | 6107051   | 5679892   | 5265806        | 980464        | 886540           | 4379063       | 28810941 (83.9%)      | 27596750 (80.4%)  |
chr14        | 41460269      | 9678269   | 794980     | 1080314   | 7802975   | 7226647   | 6668634        | 1306986       | 1177106          | 5491258       | 33542710 (80.9%)      | 32475908 (78.3%)  |
chr15        | 41483448      | 9466026   | 767928     | 1040182   | 7657916   | 7101722   | 6562680        | 1255801       | 1133216          | 5429039       | 33585637 (81.0%)      | 32747719 (78.9%)  |
chr16        | 36990834      | 8523812   | 683735     | 944310    | 6895767   | 6384354   | 5889245        | 1154066       | 1038891          | 4849963       | 29820666 (80.6%)      | 28140258 (76.1%)  |
chr17        | 40844869      | 9509196   | 790475     | 1065255   | 7653466   | 7089503   | 6543522        | 1291396       | 1163555          | 5379708       | 33351569 (81.7%)      | 32297540 (79.1%)  |
chr18        | 37040606      | 8514807   | 689451     | 934855    | 6890501   | 6377321   | 5880972        | 1114153       | 1002056          | 4878380       | 29742506 (80.3%)      | 28106975 (75.9%)  |
chr19        | 31504643      | 7505226   | 606060     | 830580    | 6068586   | 5603100   | 5152963        | 1020938       | 916595           | 4236048       | 25755851 (81.8%)      | 24193841 (76.8%)  |
chr20        | 32537810      | 7308444   | 614238     | 818789    | 5875417   | 5447570   | 5033913        | 969006        | 872105           | 4161355       | 26164251 (80.4%)      | 24941362 (76.7%)  |
chr22        | 35317527      | 8501749   | 707735     | 951427    | 6842587   | 6307831   | 5791610        | 1178643       | 1054610          | 4735743       | 27880119 (78.9%)      | 26337572 (74.6%)  |
chr23        | 45375999      | 10112574  | 792755     | 1096778   | 8223041   | 7624547   | 7045510        | 1311479       | 1181745          | 5863276       | 35441094 (78.1%)      | 33272875 (73.3%)  |
all_main_chr | 876670382     | 201347802 | 16451413   | 22301042  | 162595347 | 150547040 | 138886900      | 26777895      | 24096738         | 114778666     | 700196486 (79.9%)     | 661988169 (75.5%) |
U_scaffolds  | 29421111      | 5649756   | 381473     | 519431    | 4748852   | 4347073   | 3961264        | 63695         | 56870            | 3901640       | 10757890 (36.6%)      | 2960951 (10.1%)   |
chrM         | 16797         | 5766      | 122        | 89        | 5555      | 4856      | 4219           | 73            | 49               | 4169          | 16380 (97.5%)         | 14532 (86.5%)     |

### Accessability masks

Bed format files containing accessible regions of the genome based on quality in the 255-individual QC subset can be found here: 
/rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/malawi_callset/mask/bed

Note that these *_.pass bed files do not consider HighDP/LowDP (coverage cutoffs based on all 2198 callset individuals).

### Updates

**22.11.2021.**
New bcftools-generated callset (version 1.3) available. This callset includes unplaced scaffolds and the mitochondrial genome.

**13.05.2020.**
Ancestral allele inference based on parsimony using the outgroups *Simochromis diagramma*, *Neolamprologus brichardi*, and *Oreochromis niloticus*. Updated all-variant and biallelic VCFs to version 1.2.
