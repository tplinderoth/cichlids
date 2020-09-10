#!/usr/bin/perl

# trioEval.pl

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my $vcf;
my @parents;
my @offspring;
my $mingq = 0;
my $filter;

die(qq/
trioEval.pl [arguments]
OR
cat vcf | trioEval.pl [arguments]

Arguments:
--vcf         VCF file, can be (b)gzipped
--parents     space-delimited VCF IDs of both parents
--offspring   space-delimited VCF IDs of offspring
--minGQ       Minimum genotype quality required for both parents to analyze site [$mingq]
--filter      Discard sites with any FILTER flag present in this ','-delimited list

*Writes all sites with genotyping errors to STDOUT and a summary of genotyping error
 counts to STDERR.
\n/) if (-t STDIN && !@ARGV);

# parse arguments

GetOptions('vcf=s' => \$vcf, 'parents=s{2}' => \@parents, 'offspring=s{1,}' => \@offspring, 'minGQ=f' => \$mingq, 'filter=s' => \$filter);

my $vcffh;
if ($vcf) {
	open($vcffh, '<', $vcf) or die("Unable to open VCF file $vcf: $!\n");
	sysread $vcffh, my $magic, 2;
	close $vcffh;

	if ($magic eq "\x1f\x9d") {
		open($vcffh, '<', $vcf);
	} elsif ($magic eq "\x1f\x8b") {
		$vcffh = new IO::Zlib;
		$vcffh->open($vcf, "rb");
	} else {
		die("Unrecognized compression type for $vcf\n");
	}
} else {
	die("No input VCF\n") if (-t STDIN);
	$vcffh = \*STDIN;
}

my %discard;
map {$discard{$_} = 1} split(',',$filter) if $filter;

die("--mingq must be >= 0\n") if ($mingq < 0);

# get subset of individuals to analyze ready

my %subset;
map {$subset{$_} = 1} (@parents, @offspring);

my $vcfline;
do {chomp($vcfline = <$vcffh>)} until ($vcfline =~ /^#CHROM/);
my @vcftok = split(/\s+/, $vcfline);

for (my $i=9; $i <= $#vcftok; $i++) {
	$subset{$vcftok[$i]} = $i if exists $subset{$vcftok[$i]};
}

foreach (keys %subset) { die("ERROR: $_ not found in VCF\n") if ($subset{$_} < 9); }

my @famidx; # array index of individuals in VCF starting with the two parents
map {push @famidx, $subset{$_}} (@parents, @offspring);

my %errcounts; # keys are individual's vcf array index
foreach (@offspring) {
	$errcounts{$subset{$_}}{id} = $_; # offspring name
	$errcounts{$subset{$_}}{total} = 0; # total number called genotypes
	$errcounts{$subset{$_}}{error} = 0; # total number wrong genotypes
}

# analyze VCF
$"="\t";
print STDOUT "CHR\tPOS\tREF\tALT\t@parents\t@offspring\n";

while (<$vcffh>) {
	@vcftok = split(/\s+/, $_);

	if ($filter) {
		my $filter_fail = 0;
		foreach my $f (split(';',$vcftok[6])) {
			if ($discard{$f}) {
				$filter_fail = 1;
				last;
			}
		}
		next if ($filter_fail);
	}

	# get subfields of needed info
	my $gqidx;
	my @format = split(':', $vcftok[8]);
	for ($gqidx = 0; $gqidx <= $#format; $gqidx++) {
		last if ($format[$gqidx] eq 'GQ')
	}
	die("ERROR: No genotype quality information for $vcftok[0] $vcftok[1]\n") if ($mingq && ($gqidx < 1 || $gqidx > $#format));

	# parse parental genotypes
	my @p1 = split(':', $vcftok[$famidx[0]]);
	my @p1geno = split(/\/|\|/, $p1[0]);
	my @p2 = split(':', $vcftok[$famidx[1]]);
	my @p2geno = split(/\/|\|/, $p2[0]);

	# check if parental genotypes are appropriate for analysis
	my ($p1allele, $p2allele);
	next if ($mingq && ($p1[$gqidx] eq '.' || $p2[$gqidx] eq '.' || $p1[$gqidx] < $mingq || $p2[$gqidx] < $mingq));
	next if ($p1[0] =~ /\./ || $p2[0] =~ /\./);
	
	if ($p1geno[0] == $p1geno[1]) {
		$p1allele = $p1geno[0];
	} else {
		next;
	}

	if ($p2geno[0] == $p2geno[1]) {
		$p2allele = $p2geno[0];
	} else {
		next;
	}

	# check for offspring genotype errors
	my $error = 0;
	foreach my $offidx (@famidx[2..$#famidx]) {
		my @offgeno = split(/\/|\|/, $1) if ($vcftok[$offidx] =~ /^([^:]+)/);
		next if ($offgeno[0] =~ /\./); # missing offspring genotype
		unless (($offgeno[0] == $p1allele && $offgeno[1] == $p2allele) || ($offgeno[0] == $p2allele && $offgeno[1] == $p1allele)) {
			$errcounts{$offidx}{error}++;
			$error = 1;
		}
		$errcounts{$offidx}{total}++;
	}

	print STDOUT "@vcftok[0,1,3,4,@famidx]\n" if ($error);
}

close $vcffh;

print STDERR "sample\tgenotypes\terrors\n";
foreach (@famidx[2..$#famidx]) {
	print STDERR "$errcounts{$_}{id}\t$errcounts{$_}{total}\t$errcounts{$_}{error}\n";
}

exit;
