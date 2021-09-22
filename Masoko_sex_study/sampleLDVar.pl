#!/usr/bin/perl

# sampleLDVar.pl

use warnings;
use strict;
use Getopt::Std;

my ($snpfile, $plinkprefix, $outprefix);
my $minld = 0.5;
my $maxdist = 1000000;
my $focal = undef;
my $nsamp = 0;

die(qq/
sampleLDVar [options] <snpfile> <plinkprefix> <outprefix>

Input:
snpfile: 1-column file of SNPs (1 per row) <STRING>
plinkprefix: Prefix of plink *.bed, *.bim, *.fam files <STRING>
outprefix: Output file prefix <STRING>

Options:
-a <FLOAT>: Minimum r^2 with focal SNP [$minld]
-b <INT>: Maximum distance in base pairs (bp) from focal SNP to calculate r^2 [$maxdist]
-c <STRING>: comma-delimited list of focal SNP IDs
-n <INT>: Number of windows to randomly sample [All SNPs in snpfile]

Description:
Calculates the standard deviation in the distance between a focal SNP with all other SNPs
in a window having r^2 greater than a cutoff (-a).
\n/) if (!@ARGV || scalar(@ARGV) < 3);

$outprefix = pop @ARGV;
$plinkprefix = pop @ARGV;
$snpfile = pop @ARGV;

my %opts = (a => undef, b => undef);
getopts('a:b:c:n:', \%opts);
$minld = $opts{a} if $opts{a};
$maxdist = $opts{b} if $opts{b};
my @focalarr = split(',',$opts{c}) if $opts{c};
$nsamp = $opts{n} if $opts{n};

open(my $snpfh, '<', $snpfile) or die("Unable to open snp position file $snpfile: $!\n");
die("ERROR: Unable to find plink bedfile ${plinkprefix}.bed") if (!-e "$plinkprefix.bed");
my $outname = "${outprefix}.ldvar";
open(my $outfh, '>', $outname) or die("Unable to open output file $outname: $!\n");
print $outfh "CHR\tPOS\tN_SNPS\tMEAN_DIST\tSTDV_DIST\n";

my $nsamptext = $nsamp > 0 ? $nsamp : 'all';
print STDERR "\nsnp position file: $snpfile\nplink files prefix: $plinkprefix\noutput file: $outname\nminimum r2: $minld\nmaximum distance: $maxdist\nNumber SNPs to sample: $nsamptext\n";
print STDERR "Focal SNP: ",  join(',',@focalarr), "\n" if (@focalarr);
print STDERR "\n";

my %chrlen = (
chr1 => 41162407,
chr2 => 38215696,
chr3 => 52512415,
chr4 => 33433079,
chr5 => 38678279,
chr6 => 41428020,
chr7 => 68397778,
chr8 => 25827708,
chr9 => 35913139,
chr10 => 34074121,
chr11 => 36676067,
chr12 => 38669361,
chr13 => 33145951,
chr14 => 39985354,
chr15 => 40222765,
chr16 => 35801853,
chr17 => 39448915,
chr18 => 35852588,
chr19 => 30863130,
chr20 => 31467755,
chr22 => 35011464,
chr23 => 43921804);

my %exclude;
my @sites;
if (@focalarr) {
	foreach(@focalarr) {
		$exclude{$_} = 1;
	}
	@sites = @focalarr;
}

while(<$snpfh>) {
	chomp;
	push @sites, $_ unless (exists $exclude{$_});
}
close $snpfh;

my $ldkb = $maxdist/1000;
my $i = 0 - scalar(@focalarr);
my $n_na = 0;

while ($i < $nsamp) {
	die("Fewer remaining sites in snp list than snps to test\n") if (scalar(@sites) == 0);
	my $draw;
	if ($i < 0) {
		$draw = 0;
	} else {
		$draw = int(rand($#sites));
	}
	my $snp = $sites[$draw];
	my @pos = split(':',$snp);
	if ($pos[1]-$maxdist > 0 && $pos[1] + $maxdist <= $chrlen{$pos[0]}) {
		my $rv = system("plink --bfile $plinkprefix --r2 --ld-snp $snp --ld-window-kb $ldkb --ld-window $maxdist --ld-window-r2 0 --out $outprefix\n");
		if ($rv != -1) {
			open(my $plinkfh, '<', "${outprefix}.ld") or die("Unable to open plink LD file ${outprefix}.ld: $!\n");
			<$plinkfh>; # skip header
			# use Welford's online algorithm for calculating the distance variance
			my ($n, $mean, $msq, $mean_prev, $s) = (0, 0, 0, 0, undef);
			while (my $plinkline = <$plinkfh>) {
				chomp($plinkline);
				my @l = split(/\s+/, $plinkline);
				if ($l[7] > $minld) {
					my $dist = abs($l[2]-$l[5]);
					$n++;
					$mean_prev = $mean;
					$mean += ($dist - $mean)/$n;
					$msq += ($dist - $mean)*($dist - $mean_prev);
				}
				if ($n > 1) {
					$s = sqrt($msq/($n-1)); # sample standard deviation
					# $s = sqrt($msq/$n); # population standard deviation
				}
			}
			if ($s) {
				print $outfh "$pos[0]\t$pos[1]\t$n\t$mean\t$s\n";
				$i++;
			} else {
				$n_na++;
				die("ERROR: Unable to calculate variance for focal SNP $snp\n") if ($i < 0);
			}
			unlink "${outprefix}.ld";
			unlink "${outprefix}.log";
			unlink "${outprefix}.nosex";
		} else {
			die "Failure to run plink: $!\n";
		}
	} else {
		print STDERR "Skipping $snp due to edge effects\n";
	}
	splice(@sites, $draw, 1);
}

close $outfh;

print STDERR "\nFinished\nNumber SNPs with NA variance: $n_na\n";

exit;
