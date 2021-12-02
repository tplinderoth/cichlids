#!/usr/bin/perl

# genoCheck.pl

use warnings;
use strict;

die(qq/
genoCheck.pl <vcf>
OR
cat <vcf> | genoCheck.pl

Quickly checks whether every line contains the expected number of individual entries.
\n/) if (!@ARGV && -t STDIN);

my $vcf = pop @ARGV if (@ARGV);
my $vcffh;
if ($vcf) {
	open($vcffh, '<', $vcf) or die("Unable to open VCF file $vcf: $!\n");
	sysread $vcffh, my $magic, 2;
	close $vcffh;

	if ($magic eq "\x1f\x8b") {
		# gzipped file
		$vcffh = new IO::Zlib;
		$vcffh->open($vcf, "rb");
	} else {
		open($vcffh, '<', $vcf);
	}
} else {
	die("No input VCF\n") if (-t STDIN);
	$vcffh = \*STDIN;
}

my $line;
do {
	chomp($line = <$vcffh>);
} until (eof($vcffh) || $line !~ /^##/);

die("No field descriptor line found in VCF\n") if $line !~ /^#CHROM/;
my @tok = split(/\s+/, $line);
my $nsamples = scalar(@tok) - 9;

my @fail;
my $nfail = 0;
my $nsites = 0;
while (<$vcffh>) {
	chomp;
	my $n = scalar(split(/\s+/, $_)) - 9;
	if ($n != $nsamples) {
		push @fail, $_;
		$nfail++;
	}
	$nsites++;
}
close $vcffh;

print STDOUT "===== SUMMARY =====\n";
print STDOUT "Number sites in VCF: $nsites\n";
print STDOUT "Expected number samples in VCF: $nsamples\n";
print STDOUT "Number of sites without the expected number of individuals:  $nfail\n";

if ($nfail > 0) {
	print STDOUT "\n===== FAILED SITES =====\n";
	foreach (@fail) {
		print "$_\n";
	}
}

exit;

