#!/usr/bin/perl

# maskStarGeno.pl

# masked individuals have the following genotype information:
# GT set to missing, "./.", and all other FORMAT subfields to zero
# each PL value becomes 0, e.g. "0,0,0"
# encodes GQ for masked individuals as 0

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my $version='0.0.1';

die(qq/
maskStarGeno.pl v $version

maskStarGeno.pl <input vcf>
OR
cat VCF | maskStarGeno.pl [options]
\n/) if (!@ARGV && -t STDIN);

my $vcf = pop @ARGV;

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
my $command = "##maskStarGenoCommand=<ID=maskStarGeno.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
$command .= "$vcf" if ($vcf);
$command .= "\">";

my $line;
while (($line = <$vcffh>) =~ /^#/) {
	if ($line =~ /^#CHROM/) {
		print "$command\n";
		last;	
	}
	print $line;
}
print $line;

# process VCF sites

$" = "\t";

my $nmask = 0; # number of masked individuals
my $nsites = 0; # number of sites

while ($line = <$vcffh>) {
	my @tok = split(/\s+/, $line);
	die("Multiple reference alleles found at $tok[0] $tok[1]\n") if ($tok[3] =~ /,/);
	die("Unexpecte allele in reference at $tok[0] $tok[1]\n") if ($tok[4] =~ /[^ACGT]i/);
	$nsites++;

	if ($tok[4] =~ /\*/) {
	# position contains star allele
		my @alt = split(/,/,$tok[4]);
		my $nalleles = scalar(@alt) + 1;
		my $staridx = 1;
		foreach my $allele (@alt) {
			last if ($allele eq '*');
			$staridx++;
		}

		my $npl = ($nalleles*($nalleles+1)/2); # number of genotype likelihoods for nalleles alleles
		my $mask = "./.:" . "0," x ($nalleles-1) . "0:0:0:" . "0," x ($npl-1) . "0"; #GT:AD:DP:GQ:PL, GQ=0 indicates masked individual

		for (my $i = 9; $i <= $#tok; $i++) {
			my ($a1, $a2) = ($1, $2) if ($tok[$i] =~ /^([.|\d]+)\/([.|\d]+):/);
			if (defined $a1 && defined $a2) {
				if ($a1 eq '.' && $a2 eq '.') {
					next;
				}
				elsif ($a1 == $staridx || $a2 == $staridx) {
					$tok[$i] = $mask;
					$nmask++;
				}
			} else {
				die("ERROR: Unrecognized genotype format at $tok[0] $tok[1]: $tok[$i]\n");
			}
		}

		#$tok[4] =~ s/,\*//;
		#$tok[4] =~ s/\*,//;
	}

	$tok[5] = '.'; # remove QUAL values
	$tok[6] = '.' ; # remove FILTER annotations

	print "@tok\n";
}

print STDERR "Number of sites: $nsites\nTotal number masks: $nmask\n";

close $vcffh;

exit;
