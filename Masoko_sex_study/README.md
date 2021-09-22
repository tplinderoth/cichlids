Resources for the study "differential use of multiple genetic sex determination systems 
in divergent ecomorphs of an African crater lake cichlid" by Munby *et al*. (2021)
=======================================================================================

Commands and code used for analyses are contained in [Munby_etal_2021_code.txt](./Munby_etal_2021_code.txt)
<br>

All of the source data needed to regenerate figures and results using the R code in 
[Munby_etal_2021_code.txt](./Munby_etal_2021_code.txt) is available at <https://doi.org/10.5061/dryad.3tx95x6h0>

### Scripts used to perform analyses

[sampleLDVar.pl](./sampleLDVar.pl): Used to generate pairwise LD distributions discussed in Results section
*Antagonism between Y alleles and Admixture* and shown in Figure S2.
<br>

	sampleLDVar [options] <snpfile> <plinkprefix> <outprefix>
	
	Input:
	snpfile: 1-column file of SNPs (1 per row) <STRING>
	plinkprefix: Prefix of plink *.bed, *.bim, *.fam files <STRING>
	outprefix: Output file prefix <STRING>
	
	Options:
	-a <FLOAT>: Minimum r^2 with focal SNP [0.5]
	-b <INT>: Maximum distance in base pairs (bp) from focal SNP to calculate r^2 [1000000]
	-c <STRING>: comma-delimited list of focal SNP IDs
	-n <INT>: Number of windows to randomly sample [All SNPs in snpfile]
	
	Description:
	Calculates the standard deviation in the distance between a focal SNP with all other SNPs
	in a window having r^2 greater than a cutoff (-a).

### Contact
Tyler Linderoth, tylerp.linderoth@gmail.com
