#!/usr/bin/perl

# speciesFst.pl

use warnings;
use strict;
use Getopt::Long;

my $version = '1.2.0';

my $spfile = undef;
my $idfile = undef;
my $vcf = undef;
my $posfile = undef;
my $passonly = undef;
my $out = undef;
my $comptype = 'all';

die(qq/
speciesFst.pl v $version

speciesFst.pl [input]

Input:
--spfile    2-column file with each row specifying (1) genus, and (2) species (required)
--idfile    Tab-delimited metadata file that must have columns (1) Sanger ID, (2) genus, and (3) species in that order (required)
--vcf       VCF file with sites to analyze (required)
--passonly  Use only PASS sites from FILTER field
--posfile   Position file with (1) chr, and (2) position of sites to analyze
--out       Output file name prefix (required)
--comptype  string specifying what species comparisons to make [$comptype]

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
\n/) unless (@ARGV && scalar @ARGV >= 8);

# check whether vcftools can be found
die("ERROR: Could not execute vcftools --> check that executable is in PATH") if system("vcftools > /dev/null 2>&1") < 0;

# parse arguments
GetOptions('spfile=s' => \$spfile, 'idfile=s' => \$idfile, 'vcf=s' => \$vcf, 'passonly' => \$passonly, 
'posfile=s' => \$posfile, 'out=s' => \$out, 'comptype=s' => \$comptype);

# check that arguments are valid
open(my $spfh, '<', $spfile) or die("Could not open species file $spfile: $!\n");
open(my $idfh, '<', $idfile) or die("Could not open species table file $idfile: $!\n");

my $vcffh;
open($vcffh, '<', $vcf) or die("Couldn't open position file $vcf: $!\n");
sysread $vcffh, my $magic, 2;
close $vcffh;
my $gz = $magic eq "\x1f\x8b" ? 1 : 0;

if ($posfile) {
	die("Could not find position file $posfile: $!\n") unless (-e $posfile);
}

my $outname = "${out}.fst";
open(my $outfh, '>', $outname) or die("Couldn't open output file $outname: $!\n");
close $outfh;

my $out_avg = $out;
$out_avg .= '.avgfst';
open(my $outavg_fh, '>>', $out_avg) or die("Couldn't open genome-wide average Fst output $out_avg: $!\n");

my $out_log;
$out_log = "${out}.log";

# read in species
my @sp;

while (<$spfh>) {
	chomp;
	$_ =~ s/\s+/_/g;
	$_ =~ s/\.+|"+|'+//g;
	push @sp, lc($_);
}

close $spfh;

die("Fst calculation requires at least 2 populations\n") if scalar(@sp) < 2;

# read in species table metadata

my %id;

<$idfh>; # skip header

while(<$idfh>) {
	chomp;
	my @tok = split('\t', $_);
	my $sp = "$tok[1]_$tok[2]";
	$sp =~ s/\s+/_/g;
	$sp =~ s/\.+|"+|'+//g;
	push @{$id{lc($sp)}}, $tok[0];
}

close $idfh;

# decide what kind of comparison to do

my (@idx1, @idx2);
my $pairwise = 1;

if ($comptype =~ /multi/) {
	$pairwise = 0;
	if ($comptype eq 'multi') {
		@idx1 = 0 .. $#sp;
	} elsif ($comptype =~ /\[([\d+,?]+)\]/) {
		@idx1 = split(/,/, $1);
		pop @idx1 if ($idx1[$#idx1] eq ',');
		@idx1 = map {$_-1} @idx1;
	} else {
		die("Unrecognized --comptype argument $comptype\n");
	}
} elsif ($comptype !~ /,/ && $comptype eq 'all') {
	@idx1 = 0 .. $#sp;
	@idx2 = 1 .. $#sp;
} elsif ($comptype =~ /(\d+),down/) {
	$idx1[0] = $1-1;
	@idx2 = $1 .. $#sp;
} elsif ($comptype =~ /^(all|\d+|\[[\d+,?]+\]),(all|\d+|\[[\d+,?]+\])$/) {
	my $first = $1;
	my $second = $2;

	if ($first eq 'all') {
		@idx1 = 0 .. $#sp;
	} elsif ($first =~ /^(\d+)/) {
		$idx1[0] = $1-1;
	} elsif ($first =~ /^\[([\d+,?]+)\]/) {
		@idx1 = split(/,/, $1);
		pop @idx1 if ($idx1[$#idx1] eq ',');
		@idx1 = map {$_-1} @idx1;
	} else {
		die("Unrecognized --comptype argument $comptype\n");
	}

	if ($second eq 'all') {
		@idx2 = 0 .. $#sp;
	} elsif ($second =~ /^(\d+)/) {
		$idx2[0] = $1-1;
	} elsif ($second =~ /^\[([\d+,?]+)\]/) {
		@idx2 = split(/,/, $1);
		pop @idx2 if ($idx2[$#idx2] eq ',');
		@idx2 = map {$_-1} @idx2;
	} else {
		die("Unrecognized --comptype argument $comptype\n");
	}
	
} else {
	die("Unrecognized --comptype argument $comptype\n");
}

if (!$pairwise) {
# multiple species Fst
	print $outavg_fh "mean_fst\tweighted_fst\n";

	my $vcftools_command = 'vcftools ';
	$vcftools_command .= $gz ? "--gzvcf $vcf " : "--vcf $vcf ";
	$vcftools_command .= "--out $out ";
	$vcftools_command .= "--remove-filtered-all " if ($passonly);
	$vcftools_command .= "--positions $posfile " if ($posfile);

	my @listfiles;

	print STDERR "calculating Fst among:\n";
	foreach my $i (@idx1) {
		print STDERR "$sp[$i]\n";		

		# make sample list files
		my $splist = "${out}_${sp[$i]}_tmp.txt";
		open(my $spfh, '>', $splist) or die("Couldn't open temporary sample file list $splist: $!\n");
		foreach my $sangerid (@{$id{$sp[$i]}}) {
			print $spfh "$sangerid\n";
		}
		close $spfh;

		$vcftools_command .= "--weir-fst-pop $splist ";
		push @listfiles, $splist;
	}

	$vcftools_command .= "2>>$out_log";

	# calulate Fst
	my $rv = system($vcftools_command);
	die("ERROR: Failure running vcftools\n") if $rv;
	rename "${out}.weir.fst", "${out}.fst";

	# extract genome-wide average Fst estimates from vcftools log file
	chomp(my $unweighted = `tail -n 4 $out_log | grep "mean Fst estimate" | cut -f 7 -d ' '`);
	chomp(my $weighted = `tail -n 4 $out_log | grep "weighted Fst estimate" | cut -f 7 -d ' '`);
	print $outavg_fh "$unweighted\t$weighted\n";

	# delete temporary files
	foreach my $tmpfile (@listfiles) {
		unlink($tmpfile) or warn("Couldn't delete temporary file $tmpfile");
	}

} else {
# process species pairs for pairwise Fst

	print $outavg_fh "species1\tspecies2\tmean_fst\tweighted_fst\n";

	my %compared;
	my $postmp = "";
	my $fstfiles;

	foreach my $i (@idx1) {

		# make sample list 1
		my $splist1 = "${out}_sp1_${sp[$i]}_tmp.txt";
		open(my $spfh, '>', $splist1) or die("Couldn't open temporary sample file list $splist1: $!\n");
		foreach my $sangerid (@{$id{$sp[$i]}}) {
			print $spfh "$sangerid\n";
		}
		close $spfh;

		for my $j (@idx2) {

			if ($j == $i || exists $compared{"${i}_${j}"} || exists $compared{"${j}_${i}"}) {
				next;
			}

			# make sample list 2
			my $splist2 = "${out}_sp2_${sp[$j]}_tmp.txt";
			open(my $spfh, '>', $splist2) or die("Couldn't open temporary sample file list $splist2: $!\n");
			foreach my $sangerid (@{$id{$sp[$j]}}) {
				print $spfh "$sangerid\n";
			}
			close $spfh;

			# calculating pairwise Fst
		
			my $pair = "${sp[$i]}_v_${sp[$j]}";
			my $fstout = "${out}_${pair}_tmp";
			my $fstval = "${fstout}.weir.fst";
		
			my $vcftools_command = 'vcftools ';
			$vcftools_command .= $gz ? "--gzvcf $vcf " : "--vcf $vcf ";
			$vcftools_command .= "--out $fstout ";
			$vcftools_command .= "--remove-filtered-all " if ($passonly);
			$vcftools_command .= "--positions $posfile " if ($posfile);
			$vcftools_command .= "--weir-fst-pop $splist1 --weir-fst-pop $splist2";
			$vcftools_command .= " 2>>$out_log";

			print STDERR "$sp[$i] vs $sp[$j]\n";
			my $rv = system($vcftools_command);
			die("ERROR: Failure running vcftools\n") if ($rv);

			if (!$postmp) {
			# extract chromosome and positions
				$postmp = $fstout . '.pos';
				$rv = system("cut -f1,2 $fstval > $postmp");
				die("ERROR: Could not extract positions from fst file $fstval\n") if ($rv);
				$fstfiles = $postmp;
			}

			# extract fst values and add header
			my $fstval2 = $fstval . 'only';
			$rv = system("{ echo '$pair' ; cut -f3 $fstval | tail -n +2; } > $fstval2");
			die("ERROR: Could not extract Fst values from file $fstval\n") if ($rv);
			$fstfiles .= " $fstval2";
			unlink($fstval) or warn("Could not delete temporary fst file $fstval");

			# extract genome-wide average Fst estimates from vcftools log file
			chomp(my $unweighted = `tail -n 4 $out_log | grep "mean Fst estimate" | cut -f 7 -d ' '`);
			chomp(my $weighted = `tail -n 4 $out_log | grep "weighted Fst estimate" | cut -f 7 -d ' '`);
			print $outavg_fh "$sp[$i]\t$sp[$j]\t$unweighted\t$weighted\n";

			unlink($splist2) or warn("WARNING: Couldn't remove temporary file $splist2");

			$compared{"${i}_${j}"} = 1;
		}
	
	unlink($splist1) or warn("WARNING: Couldn't remove temporary file $splist1");
	}

	# combine pairwise Fst results into one file
	my $combrv = system("paste $fstfiles > $outname");
	die("ERROR: Failed to combine pairwise Fst files\n") if $combrv;

	foreach my $tmpfile (split(/\s/, $fstfiles)) {
		unlink($tmpfile) or warn("Could not delete temporary file $tmpfile");
	}
}

close $outavg_fh;

exit;
