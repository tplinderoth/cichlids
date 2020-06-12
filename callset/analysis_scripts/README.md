Scripts for analyzing the cichlid callset
=========================================

## callsetTools.pl

callsetTools.pl is a collection of utilities for querying, manipulating, and analyzing VCFs. It's current functionality includes
* Obtaining Sanger IDs for samples belonging to a group such as genus, species, location, etc.
* Calculating allele frequencies (+ REF/ALT/ANC alleles, quality, and missing data info) for any group or list of groups, e.g. species, at single or multiple sites.
* Converting VCF (includines subsets) to other formats, e.g. BEAGLE, ANGSD geno, etc.

### Example usage

To obtain information for the main functions run without arguments or --help.

	./analysis_scripts/callsetTools.pl
	
	callsetTools.pl 1.2.0
	
	Requires bcftools and vcftools to be installed and in $PATH
	
	callsetTools.pl [command]
	
	commands:
	sampleID     retrieve Sanger IDs for samples belonging to a group
	alleleFreq   calculate allele frequencies
	multiFreq    calculate allele frequencies for a list of sites
	convert      convert VCF to other formats (allows subsetting)
	
To obtain information on how to run a specific command, run without arguments.

	./analysis_scripts/callsetTools.pl alleleFreq
	
	callsetTools.pl alleleFreq [chr number] [position] [sample list] [vcf]
	
	Sample list:
	2-column file in sampleID output format where each row lists (1) Sanger ID of sample, (2) group identity.
	Assumes a header line is present. Allele frequencies will be calculated for each group.
	
	VCF:
	Indexed VCF file.

### Example for calculating allele frequencies for a specific set of species.
First, gather the Sanger IDs for individuals within each species that you want to include in the allele frequency calculation:

	./callsetTools.pl sampleID species ./example_files/species_list_fst.txt ./example_files/malawi_cichlids_v1_sample_table.txt > ./example_files/id_list_fst.txt


Now calculate allele frequencies for each species:

	VCF='~/rds/rds-rd109-durbin-group/projects/cichlid/cichlid_g.vcf/main_set_2019-07-10/fAstCal1.2/GenotypeCorrected/depth_filter/malawi_variants/biallelic1.2/malawi_cichlid_v1.2_biallelic_chr23.vcf.gz'
	
	./callsetTools.pl alleleFreq 23 7998908 ./example_files/id_list_fst.txt $VCF > ./example_files/kita_chr23_7998908.freq

	cat ./example_files/kita_chr23_7998908.freq
	chr23	7998908	ref=T	alt=G	anc=T	PASS	AC=1557;AF=0.354186;AN=4396;DP=34607;NS=2198;MAF=0.354186;ExcHet=1;AA=T;VT=snp
	group	alt_freq	n	ndata
	Astatotilapia_calliptera	0	733	732
	Alticorpus_peterdaviesi	0.239130434782609	23	23
	Aulonocara_blue_orange	0	11	11
	Aulonocara_slender_orange_back	0.333333333333333	12	12
	Aulonocara_yellow	0.0769230769230769	13	13
	Aulonocara_stuartgranti	0	10	10
	Chilotilapia_rhoadesii	0.0277777777777778	18	18
	Copadichromis_chrysonotus	0.422727272727273	110	110
	Copadichromis_virginalis	0.97463768115942	138	138
	Cynotilapia_afra	0.895522388059702	67	67
	Dimidiochromis_kiwinge	0.0454545454545455	11	11
	Diplotaxodon_macrops_black_dorsal	0.34375	16	16
	Diplotaxodon_limnothrissa	0.757575757575758	33	33
	Fossorochromis_rostratus	0.223684210526316	38	38
	Labeotropheus_fuelleborni	0.611111111111111	27	27
	Labeotropheus_trewavasae	0.810344827586207	30	29
	Lethrinops_gossei	0.025	20	20
	Lethrinops_lethrinus	0.0357142857142857	14	14
	Maylandia_pearly	0	17	17
	Maylandia_callainos	0	25	25
	Maylandia_emmiltos	0.416666666666667	30	30
	Maylandia_fainzilberi	0.72972972972973	74	74
	Maylandia_zebra	0.93125	80	80
	Melanochromis_auratus	1	19	19
	Mylochromis_subocularis	0.227272727272727	22	22
	Otopharynx_argyrosoma	0.0405405405405405	37	37
	Petrotilapia_genalutea	0.933333333333333	15	15
	Petrotilapia_microgalana	0.708333333333333	12	12
	Protomelas_ornatus	0.607142857142857	14	14
	Rhamphochromis_longiceps	0	18	18
	Trematocranus_placodon	0.652777777777778	36	36
	Tropheops_chilumba	1	16	16
	Tropheops_yellow_gular	1	16	16

If you want to calculate allele frequencies for all of the above species over multiple sites at once, here's an example of how:

	./callsetTools.pl multiFreq -V ./example_files/malawi_biallelic1.2_vcfs.txt ./example_files/malawi_radiation_fst_sites.pos ./example_files/id_list_fst.txt > ./example_files/malawi_radiation_fst_sites.freq
	
	cat ./example_files/malawi_radiation_fst_sites.freq
	CHR	POS	REF	ALT	AA	FILTER	INFO	Astatotilapia_calliptera	Alticorpus_peterdaviesi	Aulonocara_blue_orange	Aulonocara_slender_orange_back	Aulonocara_yellow	Aulonocara_stuartgranti	Chilotilapia_rhoadesii	Copadichromis_chrysonotus	Copadichromis_virginalis	Cynotilapia_afra	Dimidiochromis_kiwinge	Diplotaxodon_macrops_black_dorsal	Diplotaxodon_limnothrissa	Fossorochromis_rostratus	Labeotropheus_fuelleborni	Labeotropheus_trewavasae	Lethrinops_gossei	Lethrinops_lethrinus	Maylandia_pearly	Maylandia_callainos	Maylandia_emmiltos	Maylandia_fainzilberi	Maylandia_zebra	Melanochromis_auratus	Mylochromis_subocularis	Otopharynx_argyrosoma	Petrotilapia_genalutea	Petrotilapia_microgalana	Protomelas_ornatus	Rhamphochromis_longiceps	Trematocranus_placodon	Tropheops_chilumba	Tropheops_yellow_gular
	chr23	7998908	T	G	T	PASS	AC=1557;AF=0.354186;AN=4396;DP=34607;NS=2198;MAF=0.354186;ExcHet=1;AA=T;VT=snp	0	0.239130434782609	0	0.333333333333333	0.0769230769230769	0	0.0277777777777778	0.422727272727273	0.97463768115942	0.895522388059702	0.0454545454545455	0.34375	0.757575757575758	0.223684210526316	0.611111111111111	0.810344827586207	0.025	0.0357142857142857	0	0	0.416666666666667	0.72972972972973	0.93125	1	0.227272727272727	0.0405405405405405	0.933333333333333	0.708333333333333	0.607142857142857	0	0.652777777777778	1	1
	chr15	1508395	G	A	G	PASS	AC=1422;AF=0.323623;AN=4394;DP=33803;NS=2197;MAF=0.323623;ExcHet=1;AA=G;VT=snp,widel	0	0.978260869565217	1	0.916666666666667	0.230769230769231	1	1	0.0409090909090909	0	0.865671641791045	0	0	0.0151515151515152	1	1	1	1	0	1	0.875	0.383333333333333	0.655405405405405	0.9375	0.552631578947368	0.0476190476190476	0.0540540540540541	0.533333333333333	0.916666666666667	0.178571428571429	0	1	0.15625	0.1875
	chr18	33845586	A	G	G	PASS	AC=1469;AF=0.336156;AN=4370;DP=32296;NS=2185;MAF=0.336156;ExcHet=1;AA=G;VT=snp,widel	0.00278164116828929	0	0.909090909090909	0	0	0.5	1	0	0.0217391304347826	0.932835820895522	1	0	0	1	0.740740740740741	0.566666666666667	0	0.892857142857143	1	1	1	0.97972972972973	1	0.0789473684210526	0.977272727272727	0.581081081081081	0.733333333333333	0.458333333333333	0.964285714285714	0	0.930555555555556	0.0625	0.03125
	chr23	103881	G	A	G	PASS	AC=1293;AF=0.316292;AN=4088;DP=24339;NS=2044;MAF=0.316292;ExcHet=1;AA=G;VT=snp,widel	0	1	0.954545454545455	1	0.769230769230769	1	0.0277777777777778	0.259090909090909	0.989130434782609	0.0373134328358209	0	0.807692307692308	1	0.0394736842105263	0.925925925925926	0.775862068965517	1	0	0.766666666666667	0.28	0	0.00675675675675676	0.04375	0.75	0	0.0135135135135135	0.892857142857143	1	0	1	0	0.884615384615385	1
	chr10	22379578	A	G	A	PASS	AC=1632;AF=0.372773;AN=4378;DP=32834;NS=2189;MAF=0.372773;ExcHet=1;AA=A;VT=snp,widel	0	0	0.954545454545455	0	0.461538461538462	0.95	1	0.995412844036697	0.992753623188406	0.00746268656716418	1	0.9375	0.875	1	0.981481481481482	0.983333333333333	0	1	0	0	0	0.00675675675675676	0.00625	0	0.431818181818182	1	0	0.0416666666666667	1	0.0833333333333333	0.930555555555556	0.21875	0.09375
	chr16	35428019	G	A	G	PASS	AC=1421;AF=0.324429;AN=4380;DP=31092;NS=2190;MAF=0.324429;ExcHet=1;AA=G;VT=snp	0.0116918844566713	0.826086956521739	1	1	1	1	0.361111111111111	0.00458715596330275	0.0364963503649635	0.970149253731343	0.227272727272727	0	0	0.881578947368421	0.037037037037037	0	0.925	0	1	1	1	0.993150684931507	1	0	0.772727272727273	0.0540540540540541	0.333333333333333	0.208333333333333	0.857142857142857	0	0.971428571428571	0.25	0.65625
	chr17	22387553	T	A	T	PASS	AC=1532;AF=0.348657;AN=4394;DP=31640;NS=2197;MAF=0.348657;ExcHet=1;AA=T;VT=snp,widel	0	0	0	0.0416666666666667	0.0769230769230769	0	0	0.986363636363636	0.992753623188406	0	0.954545454545455	1	0.984848484848485	0	0.0185185185185185	0.0333333333333333	0	0	1	1	0.7	0.912162162162162	0.80625	0.157894736842105	0.0909090909090909	0.918918918918919	0.833333333333333	0.958333333333333	0.0357142857142857	1	0	0.09375	0.03125
	chr10	4272020	C	A	N	PASS	AC=1270;AF=0.289294;AN=4390;DP=30287;NS=2195;MAF=0.289294;ExcHet=1;AA=N;VT=snp,widel	0.00545702592087312	0.978260869565217	0.545454545454545	1	0.384615384615385	0.65	0.861111111111111	0.0909090909090909	0	0.880597014925373	0.454545454545455	0	0	1	0	0.0166666666666667	0.1	1	0.852941176470588	0.82	0.983333333333333	0.918918918918919	0.98125	0	0.476190476190476	0.918918918918919	0	0	0.142857142857143	0	0.888888888888889	0	0
	chr6	35629870	A	G	G	PASS	AC=1660;AF=0.377273;AN=4400;DP=34558;NS=2200;MAF=0.377273;ExcHet=1;AA=G;VT=snp	0	0	0.0454545454545455	0	0	0.15	0.944444444444444	0.968181818181818	0.956521739130435	0.873134328358209	0.909090909090909	0	0	0	1	1	0	0.357142857142857	1	1	1	0.97972972972973	1	0.157894736842105	0.75	0.797297297297297	0.333333333333333	0.25	0.892857142857143	0.694444444444444	0.152777777777778	0.03125	0
	chr22	1946380	T	C	T	PASS	AC=2153;AF=0.491328;AN=4382;DP=34359;NS=2191;MAF=0.491328;ExcHet=1;AA=T;VT=snp	0.394337016574586	0.260869565217391	1	1	0.961538461538462	0.5	0	0.986363636363636	0.996376811594203	0.902985074626866	0.0909090909090909	0	0	0	0.185185185185185	0.3	0.275	0	1	1	1	1	0.99375	0.868421052631579	0.5	0.959459459459459	0.2	0.0833333333333333	0	0	0	0.375	0.03125
	chr1	37988153	T	C	T	PASS	AC=1830;AF=0.436963;AN=4188;DP=28909;NS=2094;MAF=0.436963;ExcHet=1;AA=T;VT=snp	0.000685871056241427	1	0.666666666666667	1	0.807692307692308	0.75	0	0.995454545454545	1	0.985074626865672	0	0.09375	0.0606060606060606	0	0.925925925925926	0.931034482758621	1	0	1	1	1	1	0.9875	0.842105263157895	0	0.145161290322581	0.8	0.166666666666667	0	0	0.287878787878788	0.46875	0.59375
	chr6	35466151	C	T	C	PASS	AC=1035;AF=0.236517;AN=4376;DP=33588;NS=2188;MAF=0.236517;ExcHet=1;AA=C;VT=snp	0	0.978260869565217	0.272727272727273	0.333333333333333	0.846153846153846	0.3	0.972222222222222	0.127272727272727	0	0.00746268656716418	0.227272727272727	0	0	1	0.944444444444444	0.966666666666667	0.675	0.892857142857143	0	0	0	0	0.00714285714285714	0	1	0.364864864864865	0.366666666666667	0.666666666666667	1	0.0555555555555556	1	0	0.21875
	chr23	96492	C	G	G	PASS	AC=1789;AF=0.413164;AN=4330;DP=30892;NS=2165;MAF=0.413164;ExcHet=1;AA=G;VT=snp	0.262829403606103	1	0.954545454545455	1	0.653846153846154	1	0	0.257009345794392	0.978260869565217	0.0573770491803279	0	1	1	0	0.925925925925926	0.766666666666667	1	0	0.823529411764706	0.38	0	0.00735294117647059	0.069620253164557	0.894736842105263	0.0227272727272727	0.0571428571428571	0.966666666666667	1	0	1	0.0416666666666667	0.96875	1
	chr14	15371583	T	C	C	PASS	AC=1798;AF=0.408822;AN=4398;DP=31482;NS=2199;MAF=0.408822;ExcHet=1;AA=C;VT=snp,widel	0.000683060109289617	0.91304347826087	0.909090909090909	0.791666666666667	0.923076923076923	0.65	0	0.0909090909090909	0.989130434782609	1	1	1	1	0	0.981481481481482	1	0.975	0	1	1	0.866666666666667	1	0.94375	0.842105263157895	0.0909090909090909	0	0.333333333333333	0.333333333333333	0	1	0	0.0625	0.15625
	chr20	359689	C	T	T	PASS	AC=1535;AF=0.349499;AN=4392;DP=26203;NS=2196;MAF=0.349499;ExcHet=1;AA=T;VT=snp	0	0.804347826086957	0.727272727272727	0.708333333333333	0.884615384615385	0	0.25	0.00454545454545455	0.0181159420289855	1	0.0454545454545455	1	0.863636363636364	0	0.0185185185185185	0	1	0.107142857142857	1	1	0.966666666666667	0.986486486486487	1	1	0.340909090909091	0.0405405405405405	1	0.25	0.75	1	0.0571428571428571	0.9375	0.625
	chr10	21765805	G	A	G	PASS	AC=1047;AF=0.238063;AN=4398;DP=32367;NS=2199;MAF=0.238063;ExcHet=1;AA=G;VT=snp	0.00886766712141883	0.978260869565217	0.681818181818182	1	0.538461538461538	0.75	0	0	0.0036231884057971	0.947761194029851	0	0	0	0.027027027027027	1	0.983333333333333	1	0.0357142857142857	0.852941176470588	0.7	1	0.932432432432432	0.95	0.0789473684210526	0	0.0135135135135135	0.266666666666667	0.458333333333333	0	0	0	0.0625	0.03125
	chr4	238467	T	C	C	PASS	AC=1677;AF=0.38131;AN=4398;DP=33492;NS=2199;MAF=0.38131;ExcHet=1;AA=C;VT=snp,widel	0.00614754098360656	1	0.0909090909090909	1	0.961538461538462	0	0.722222222222222	0.159090909090909	0.905797101449275	0	0.181818181818182	1	1	0.105263157894737	0.037037037037037	0.0166666666666667	1	0	0.264705882352941	0.02	0.05	0.912162162162162	0.0125	0.447368421052632	0.818181818181818	0.662162162162162	1	1	0.285714285714286	1	0.763888888888889	1	1
	chr7	355601	G	A	G	PASS	AC=1253;AF=0.285162;AN=4394;DP=31162;NS=2197;MAF=0.285162;ExcHet=1;AA=G;VT=snp	0.00889192886456908	0.782608695652174	0.909090909090909	0.125	0.884615384615385	1	0	0	0.0181159420289855	0.917910447761194	0	0.375	0	0.0131578947368421	0.981481481481482	1	1	0.5	1	1	0.9	0.993243243243243	0.925	0.631578947368421	0	0.0405405405405405	0.133333333333333	0.583333333333333	0	0	0	0.9375	0.96875
	chr12	24555031	T	G	G	PASS	AC=1601;AF=0.364526;AN=4392;DP=32912;NS=2196;MAF=0.364526;ExcHet=1;AA=G;VT=snp	0.000684931506849315	0.195652173913043	0.727272727272727	0	0.0769230769230769	0.444444444444444	1	0.968181818181818	0.949275362318841	0	1	0.125	0.287878787878788	1	0	0	0.025	0.964285714285714	0	0	0	0.0743243243243243	0.00625	1	1	0.918918918918919	0.1	0.0416666666666667	1	1	1	0.875	0.71875
	chr16	8231407	G	A	G	PASS	AC=1736;AF=0.396528;AN=4378;DP=30462;NS=2189;MAF=0.396528;ExcHet=1;AA=G;VT=snp	0.000687757909215956	0.0869565217391304	0.227272727272727	0	0.846153846153846	1	1	1	0.996323529411765	0.0223880597014925	0.954545454545455	0.75	0.96969696969697	0.605263157894737	0	0.0333333333333333	0.725	1	0	0	0.0166666666666667	0.0135135135135135	0	0	1	1	0.133333333333333	0	1	1	0.972222222222222	0.65625	0.21875

## looksite.pl

looksite.pl is a tool for querying SNPs from the callset VCFs.<br/>
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

How to set the VCFDIR, METADATA, and HEADER links on the Cambridge CSD3 (assuming you are in this 'analysis_scripts' directory containing looksite.pl):

	ln -s ~/rds/rds-rd109-durbin-group/projects/cichlid/cichlid_g.vcf/main_set_2019-07-10/fAstCal1.2/GenotypeCorrected/depth_filter/malawi_variants/vcf1.2 VCFDIR

	ln -s ../200128.cichlid_callset_metadata.txt METADATA

	ln -s ../HEADER HEADER

## speciesFst.pl

speciesFst.pl is a wrapper around [VCFtools](https://vcftools.github.io/) to calculate per-site and genome-average Fst for all pairwise combinations of species provided in a list. The resulting Fst values are all nicely bundled together in two files.

	perl speciesFst.pl 
	
	speciesFst.pl v 1.2.0
	
	speciesFst.pl [input]
	
	Input:
	--spfile    2-column file with each row specifying (1) genus, and (2) species (required)
	--idfile    Tab-delimited metadata file that must have columns (1) Sanger ID, (2) genus, and (3) species in that order (required)
	--vcf       VCF file with sites to analyze (required)
	--passonly  Use only PASS sites from FILTER field
	--posfile   Position file with (1) chr, and (2) position of sites to analyze
	--out       Output file name prefix (required)
	--comptype  string specifying what species comparisons to make [all]
	
	Genus and species names in spfile should match those in the idfile.
	
	comptype arguments, where INT refers to the species in the INT row from the top of the spfile:
	multi      Fst among all species (not pairwise)
	SET,multi  Fst among all species in SET = '[INT1,INT2,...,INTn]' (without the single quotes)
	all        all pairwise comparisons
	INT,all    pairwise comparisons between the INT species and all others
	INT,down   pairwise comparisons between the INT species and all others below it
	INT1,INT2  pairwise comparison betwen INT1 and INT2 species
	
	A set of species can be specified as '[INT1,INT2,...,INTn]' (without the single quotes) in place of INT
	to specify multiple species, in which case all species in the first set will be compared to all species 
	in the second set. 
	
	Dependencies:
	VCFtools must be installed and in user's PATH
