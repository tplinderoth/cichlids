#!/usr/bin/perl

# standardizeMissingGeno.pl

# expects FORMAT = GT:AD:DP:GQ:PL

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my $version = '1.0.0';

my $noheader = undef;

die(qq/
standardizeMissingGeno.pl v$version

standardizeMissingGeno.pl <input vcf> [options]
OR
cat VCF | standardizeMissingGeno.pl [options]

options:
--noheader   Do not append command to VCF header
\n/) if (!@ARGV && -t STDIN);

GetOptions('noheader' => \$noheader);

# parse user args
my %userargs = ('--noheader' => $noheader);

my $vcf = pop @ARGV if (@ARGV);

my $vcffh;
if ($vcf) {
	open($vcffh, '<', $vcf) or die("Couldn't open VCF file $vcf: $!\n");
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

# proces VCF header

my $datestr = localtime();
my $command = "##standardizeMissingGenoCommand=<ID=standardizeMissingGeno.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
	$command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

my $line;
while ($line = <$vcffh>) {
	last if $line =~ /^#CHROM/;
	print STDOUT $line;
}
print STDOUT $command unless defined $noheader;
print STDOUT $line;

# Standardize missing genotypes

my $fmt='GT:AD:DP:GQ:PL'; # checks for this format right now but can be generalized
my $nfields = ($fmt =~ tr/://) + 1;

$" = "\t";

my $nsites = 0;
my $nstand = 0;
my $nmiss = 0;

my $gq_mask = '0';
my $dp_mask = '0';

my $gq_miss = '.';
my $dp_miss = '0';

my $chr;

while ($line = <$vcffh>) {
	my @tok = split(/\s+/, $line);
	$nsites++;
	$chr = $tok[0];

	my $refcount = 1 + ($tok[3] =~ tr/,//); # Don't expect more than 1 ref, but will count to play it safe
	my $altcount = $tok[4] eq '.' ? 0 : ($tok[4] =~ tr/,//) + 1;
	my $nalleles = $refcount + $altcount;
	my $npl = ($nalleles*($nalleles+1)/2); # number of genotype likelihoods for nalleles alleles

	my $pl_miss = "0," x ($npl-1) . "0"; # missing PL
	my $ad_miss = "0," x ($nalleles-1) . "0"; # missing AD 
	my $missing_str = "./.:$ad_miss:$dp_miss:$gq_miss:$pl_miss";
	my @miss_arr = ('./.', $ad_miss, $dp_miss, $gq_miss, $pl_miss);
	my $mask_str = "./.:$ad_miss:$dp_mask:$gq_mask:$pl_miss";

	for (my $i = 9; $i <= $#tok; $i++) {
		if ($tok[$i] =~ /^\.|\d+\/\./)
		{
		# missing genotype
			$nmiss++;
			my @genoarr = split(':', $tok[$i]);
			my $change = 0;

			if ($genoarr[0] ne './.') {
			# ensure formating of GT subfield
				$genoarr[0] = './.';
				$change = 1;
			}

			my $n_genofields = scalar(@genoarr);

			if ($n_genofields < $nfields) {
			# subfields are missing
			# vcf standard states that only trailing subfields can dropped so wil fill in accordingly
				$tok[$i] = join(':',@genoarr[0..$#genoarr]) . ":" . join(':',@miss_arr[$#genoarr+1..$#miss_arr]);
				$change = 1;
				$nstand += $change;
				next;
			}

			if ($genoarr[4] eq '.' && $genoarr[3] eq '.' && $genoarr[2] eq '.' && $genoarr[1] eq $ad_miss) {
			# adjust my old masking
				$tok[$i] = $mask_str;
				$change = 1;
				$nstand += $change;
				next;
			}

			if ($genoarr[4] eq '.') {
			# ensure formatting of PL subfield
				$genoarr[4] = $pl_miss;
				$genoarr[3] = $gq_miss;
				$change = 1;
			}

			if ($change) {
				$tok[$i] = join(':', @genoarr);
				$nstand++;
			}
		}
	}

	print "@tok\n";
}

close $vcffh;

print STDERR "last contig processed: $chr\n";
print STDERR "Number VCF sites: $nsites\n";
print STDERR "Numer uncalled genotypes: $nmiss\n";
print STDERR "Number standardized genotype fields: $nstand\n";

exit;
