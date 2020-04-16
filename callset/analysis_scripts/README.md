Scripts for analyzing the cichlid callset
=========================================

### looksite.pl

loosite.pl is a tool for querying SNPs from the callset VCFs.<br/>
It requires the HEADER file, a 'METADATA' link to the cichlid metadata table, and a 'VCFDIR' link to the VCF directory.

	perl looksite.pl
	usage: perl looksite.pl [-ref] <chr> <pos> [<field>] where field is simple_id, genus, species, clade
	    need to set HEADER, METADATA, VCFDIR

	perl looksite.pl 23 7997094
	chr23 7997094 C A AC=88;AF=0.0200273;AN=4394;DP=33744;NS=2197;MAF=0.0200273;ExcHet=1;AA=C;VT=snp
	all     88      4394

	perl looksite.pl 23 7997094 clade
	chr23 7997094 C A AC=88;AF=0.0200273;AN=4394;DP=33744;NS=2197;MAF=0.0200273;ExcHet=1;AA=C;VT=snp
	Mbuna   88      1078

	perl looksite.pl 23 7997094 genus
	chr23 7997094 C A AC=88;AF=0.0200273;AN=4394;DP=33744;NS=2197;MAF=0.0200273;ExcHet=1;AA=C;VT=snp
	Labeotropheus   1       114
	Maylandia       87      472

	perl looksite.pl 23 7997094 species
	chr23 7997094 C A AC=88;AF=0.0200273;AN=4394;DP=33744;NS=2197;MAF=0.0200273;ExcHet=1;AA=C;VT=snp
	Labeotropheus fuelleborni       1       54
	Maylandia """pearly"""  34      34
	Maylandia callainos     50      50
	Maylandia fainzilberi   2       148
	Maylandia zebra 1       160

	perl looksite.pl 23 7997094 simple_id | grep MayFai
	MayFai57        2       2

	grep MayFai57 METADATA 
	cichlid6994163  MayFai57        D06-D10 15.9    1       Mbuna           Maylandia       fainzilberi     F?      Chilumba        Chitande_island

	perl looksite.pl 23 7997094 simple_id | grep MayZeb
	MayZeb6 1       2

	grep MayZeb6 METADATA 
	CICHM16429639   MayZeb6 D01-J01 38.4    1       Mbuna           Maylandia       zebra   F       Nkhata_Bay      Viking_reef     OB

	perl looksite.pl 23 7997094 simple_id | grep LabFue
	LabFue17        1       2

	grep LabFue17 METADATA 
	cichlid6994096  LabFue17        D05-F04 16.3    1       Mbuna           Labeotropheus   fuelleborni             Chilumba        Luwino_reef

### speciesFst.pl

speciesFst.pl is a wrapper around [VCFtools](https://vcftools.github.io/) to calculate per-site and genome-average Fst for all pairwise combinations of species provided in a list. The resulting Fst values are all nicely bundled together in two files.

	perl speciesFst.pl 
	
	speciesFst.pl [input]
	
	Input:
	--spfile    2-column file with each row specifying (1) genus, and (2) species (required)
	--idfile    Tab-delimited metadata file that must have columns (1) Sanger ID, (2) genus, and (3) species in that order (required)
	--vcf       VCF file with sites to analyze (required)
	--passonly  Use only PASS sites from FILTER field
	--posfile   Position file with (1) chr, and (2) position of sites to analyze
	--out       Output file name prefix (required)
	
	Fst will be calculated between all pairs in the spfile.
	Genus and species names in spfile should match those in the idfile.
	
	Dependencies:
	VCFtools must be installed and in user's PATH
	
	
	perl speciesFst.pl --spfile species_list_fst.txt --idfile malawi_cichlids_v1_sample_table.txt --vcf malawi_cichlid_v1_biallelic_chr1.vcf.gz --passonly --out malawi_radiation_chr1
	
