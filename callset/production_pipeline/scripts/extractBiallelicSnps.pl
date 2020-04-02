#!/usr/bin/perl

# extractBiallelicSnps.pl

# notes:
# requires bcftools to be in user's path
# Assumes input VCF indels are left-aligned and normalized: bcftools norm -f <fasta>

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my $version = '1.0.0';

my $dpbounds = undef;
my $minmaf = 0;

die(qq/
extractBiallelicSnps.pl v$version

extractBiallelicSnps.pl [parameters] <input vcf>
OR
cat VCF | extractBiallelicSnps.pl [parameters]

Input parameters:
--dpbounds   LowDP and HighDP bounds to annotate FILTER filed, format '<INT lower bound>,<INT upper bound>'
--minmaf     Only keep sites with minor allele frequency greater than FLOAT [$minmaf]

Dependencies:
Requires bcftools executable to be in user PATH

Notes:
* Assumes input VCF indels are left-aligned and normalized (bcftools norm -f <fasta>)
\n/) if (!@ARGV && -t STDIN);

# check whether bcftools executable is in user's path

die("ERROR: Could not execute bcftools --> check that executable is in PATH") if system("bcftools > /dev/null 2>&1") < 0;
chomp(my @bcftools_out = qx\bcftools version\);

my $bcfv = $1 if ($bcftools_out[0] =~ /(\S+)$/);
my $htsv = $1 if ($bcftools_out[1] =~ /(\S+)$/);

my $depv = "bcftools-${bcfv}+htslib-${htsv}";

# parse user args

GetOptions('dpbounds=s' => \$dpbounds, 'minmaf=f' => \$minmaf);
my %userargs = ('--dpbounds' => $dpbounds, '--minmaf' => $minmaf);

my @dpcutoff;
if ($dpbounds) {
	@dpcutoff = split(/,/, $dpbounds);
	if (scalar(@dpcutoff) != 2) {
		die("ERROR: --dpbounds requires 2 values, a lower and an upper coverage threshold\n");
	} elsif ($dpcutoff[0] >= $dpcutoff[1]) {
		die("ERROR: --dpbounds upper bound must be greater than the lower bound value\n");
	}
}

die("--minmaf should be a value in [0,0.5]\n") if ($minmaf < 0 || $minmaf > 0.5);

my $rand = int(rand(9999));

my $vcf = pop @ARGV if (@ARGV);

my $tmpfile;
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

	$tmpfile = $vcf;
	$tmpfile =~ s/\.gz$|$/\.$rand\.tmp/;

} else {
	die("No input VCF\n") if (-t STDIN);
	$vcffh = \*STDIN;
	$tmpfile = "extractedSnps.vcf.$rand.tmp"
}

die("ERROR: Can't extract SNPs because temporary file $tmpfile already exists\n") if (-e $tmpfile);
open(my $tmpfh, '>', $tmpfile) or die("Can't open temporary output file $tmpfile: $!\n");

# declare variables to hold stats
my $nsites = 0; # number of sites in VCF
my $nsnps = 0; # total number SNPs
my $nsnp_sites; # total number sites containing SNPs
my $nsnp_widel = 0; # number of SNP sites within a deletion
my $nbisnps = 0; # number sites with biallelic snps
my $nbisnps_widel = 0; # number sites with biallelic snps within a deletion
my $ninsert = 0; # number of insertions
my $ndelete = 0; # number of deletions
my $nmask = 0; # number of times individual genotypes were masked

# proces VCF header

my $datestr = localtime();
my $command = "##extractBiallelicSnpsCommand=<ID=extractBiallelicSnps.pl,Version=${version}+$depv,Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
	$command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

# print pretty header

my @header = makeHeader($vcffh, $command, $dpbounds);

foreach my $hline (@header) {
	print $tmpfh $hline;
}

# Extract biallelic SNPs from VCF

$" = "\t";

my @indel = (0,0); # start and stop coordinates for indel
my ($end, $aa, $vartype, $chr);
while (my $line = <$vcffh>) {
	my @tok = split(/\s+/, $line);
	$nsites++;
	$chr = $tok[0];

	die("Multiple reference alleles found at $tok[0] $tok[1]") if ($tok[3] =~ /,/);
	die("Unexpecte allele in reference at $tok[0] $tok[1]") if ($tok[4] =~ /[^ACGT]i/);

	my $reflen = length($tok[3]);

	if ($reflen > 1) {
	# left-alignment means this is a deletion wrt reference
		# figure out if site is within a new deletion
		my $end = $tok[1] + ($reflen-1);
		if ($end > $indel[1]) {
			my $start = $tok[1]+1;
			$indel[0] = $tok[1]+1 if $start > $indel[1]; # this is the first base actually deleted
			$indel[1] = $end; # this is the last base actually deleted             
		}
	}
	my $widel = ($tok[1] >= $indel[0] && $tok[1] <= $indel[1]) ? 1 : 0;

	# classify the type of variant
	my @altarr = split(/,/, $tok[4]);
	my $issnp = 0;
	my $altn;
	my @snp_alleles = ();
	my $refbase = substr $tok[3], 0, 1;
	my $i = 1;

	foreach my $a (@altarr) {
		
		if ($a eq '*') {
			# do not record or keep star alleles
			$i++;
			next;
		}

		my $altlen = length($a);
		
		if ($altlen == $reflen) {
		# SNP
			$issnp = 1;
			my $altbase = substr $a, 0, 1;
			if ($altbase ne $refbase) {
				$altn = $i;
				push @snp_alleles, $altbase;
				$nsnps++;
			}
		} else {
		# indel
			$altlen > $reflen ? $ninsert++ : $ndelete++;
		}

		$i++;
	}

	$nsnp_sites++ if ($issnp);
	$nsnp_widel++ if ($issnp && $widel);

	# process biallelic SNP
	if (scalar(@snp_alleles) == 1) {
	# there is one alternate SNP allele (biallelic SNP site)
		$nbisnps++;

		# change representation of ref and kept alt allele in VCF line
		$tok[3] = $refbase;
		$altarr[$altn-1] = $snp_alleles[0];
		$tok[4] = join(',', @altarr);
		
		# mask individuals with nonref and alt snp alleles
		# could recalculate INFO subfields here while looping over genotypes but will leave that to bcftools for now
		if ($i > 2) {
			my $npl = ($i*($i+1)/2); # number of genotype likelihoods for nalleles alleles
			my $mask = "./.:" . "0," x ($i-1) . "0:0:0:" . "0," x ($npl-1) . "0"; #GT:AD:DP:GQ:PL, GQ=0 indicates masked individual
			
			for (my $j = 9; $j <= $#tok; $j++) {
				my ($a1, $a2) = ($1, $2) if ($tok[$j] =~ /^([.|\d]+)\/([.|\d]+):/);
                        
				if (defined $a1 && defined $a2) {
					if ($a1 eq '.' && $a2 eq '.') {
					next;
				}
				elsif (($a1 > 0 && $a1 != $altn) || ($a2 > 0 && $a2 != $altn)) {
					$tok[$j] = $mask;
					$nmask++;
				}
				} else {
					die("ERROR: Unrecognized genotype format at $tok[0] $tok[1]: $tok[$j]\n");
				}

			}
		}

		# adjust INFO ancestral alleles
		if ($tok[7] =~ /AA=([ACGTN*]+)/) {
			my $aa = $1;
			my $aa_base = substr $1, 0, 1;
			if ($aa_base eq $tok[3] || $aa_base eq $snp_alleles[0]) {
				$tok[7] =~ s/AA=$aa/AA=$aa_base/;
			} else {
				$tok[7] =~ s/AA=$aa/AA=N/
			}
		}

		# annotate INFO variant type
		$nbisnps_widel++ if ($widel);

		my $vt;
		if ($tok[7] =~ /VT=([^\s|;|=]+)/) {
			$vt = $1;
			my $vt_original = $vt;
			$vt .= ',snp' unless $vt =~ /snp/i;
			$vt .= ',widel' if ($widel && $vt !~ /widel/i);
			$tok[7] =~ s/VT=$vt_original/VT=$vt/;
		} else {
			$vt = 'snp';
			$vt .= ',widel' if ($widel);
			$tok[7] .= ";VT=$vt";
		}

		# remove previous dp information
		$tok[6] =~ s/PASS/\./;
		$tok[6] =~ s/LowDP//;
		$tok[6] =~ s/HighDP//;
		$tok[6] =~ s/,{2,}//; $tok[6] =~ s/^,+|,+$//g;
		$tok[6] = '.' if (!$tok[6]);

		# print site
		print $tmpfh "@tok\n";
	}

}

my $nvariants = $nsnps + $ndelete + $ninsert; # total number of variants

close $tmpfh;

# trim alt alleles not observed and recalculate INFO subfields

my $tmpfile2 = $tmpfile;
$tmpfile2 =~ s/tmp$/tmp2/;
die("ERROR: Can't extract SNPs because temporary file $tmpfile2 already exists\n") if (-e $tmpfile2);

my $bcftools_call = system("bcftools view --no-version -O u -a -U -c 1:nref $tmpfile | bcftools +fill-tags --no-version -o $tmpfile2 -- -t 'AN,AC,AF,MAF,NS,ExcHet,DP=sum(DP)'");
die("ERROR: Failure upon running bcftools\n") if ($bcftools_call);

unlink($tmpfile) or warn("WARNING: Could not remove temporary output file $tmpfile: $!\n");

# add FILTER annotations and ensure site is variable (filter on MAF)

open(my $vcffh2, '<', $tmpfile2) or die("Couldn't read temporary VCF file $tmpfile2: $!\n");

# print pretty header

@header = makeHeader($vcffh2, "", $dpbounds);
foreach my $hline (@header) {
	print STDOUT $hline;
}

my $snpvcf_lines = 0;
while (my $line = <$vcffh2>) {
	my @tok = split(/\s+/, $line);
	
	# allele frequency check that site is variable
	next if ($tok[4] eq '.');

	if ($tok[7] =~ /AF=([^\s;=]+)/) {
		my $af = $1;
		next unless ($af > 0 && $af < 1);

		# check to make sure MAF annotation is present and filter on it		
		if ($tok[7] =~ /MAF=([^\s;=]+)/) {
			my $maf = $1;
			next unless $maf > $minmaf;
		} else {
			print STDERR "WARNING: No INFO/MAF at $tok[0] $tok[1]: $tok[7]\n";
		}
	} else {
		die("ERROR: No INFO/AF at $tok[0] $tok[1]: $tok[7]\n");
	}
	
	# add DP FILTER ANNOTATIONS
	if ($dpbounds) {
		my $dp = $1 if ($tok[7] =~ /DP=([^\s;=]+)/);
		my $dptag = undef;
		if ($dp < $dpcutoff[0]) {
			$dptag = "LowDP";
		} elsif ($dp > $dpcutoff[1]) {
			$dptag = "HighDP";
		}

		if ($dptag) {
			if ($tok[6] eq "." || $tok[6] eq "PASS") {
				$tok[6] = $dptag;
			} else {
				$tok[6] .= ";${dptag}" if ($tok[6] !~ /$dptag/);
			}
		} else {
			$tok[6] = "PASS" if ($tok[6] eq ".");
		}
	}

	print STDOUT "@tok\n";
	$snpvcf_lines++;
}

unlink($tmpfile2) or warn("WARNING: Could not remove temporary output file $tmpfile2: $!\n");

# print stats

print STDERR "Last contig processed: $chr\n";
print STDERR "Number VCF sites (lines without header): $nsites\n";
print STDERR "Number variants: $nvariants\n";
print STDERR "Number SNPs: $nsnps\n";
print STDERR "Number insertions: $ninsert\n";
print STDERR "Number deletions: $ndelete\n";
print STDERR "Number SNP sites: $nsnp_sites\n";
print STDERR "Number bialleleic SNP sites: $nbisnps\n";
print STDERR "Number SNP sites within a deletion: $nsnp_widel\n";
print STDERR "Number biallelic SNP sites within a deletion: $nbisnps_widel\n";
print STDERR "Number masked genotype entries: $nmask\n";
print STDERR "Number sites in biallelic SNP VCF: $snpvcf_lines\n";

exit;

sub makeHeader {
	my ($fh, $command, $dpbounds) = @_;

	# 1000 Genomes style VT annotation: ##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
	my $vt_string = "##INFO=<ID=VT,Number=.,Type=String,Description=\"indicates what type of variant the line represents\">\n";

	my @headorder = ('fileformat', 'reference', 'contig', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'other'); # header order
	my %header = (fileformat => undef, ALT => undef, FILTER => undef, INFO => undef, FORMAT => undef, reference => undef, contig => undef, other => undef);
	my %seen = (vt => 0);
	my @headerlines;
	my $line;

	while (($line = readline($$fh)) =~ /^##/) {
		if ($line =~ /^##([^=]+)/i) {
			my $annotation = $1;

			if ($line =~ /INFO=<ID=VT,/) {
				# update VT annotation to ensure it is correct
				push @{$header{$annotation}}, $vt_string;
				$seen{vt} = 1;
			} elsif ($annotation =~ /^fileformat|ALT|FILTER|INFO|FORMAT|reference|contig/) {
				push @{$header{$annotation}}, $line;
			} else {
				push @{$header{other}}, $line;
			}
		}
	}

	push @{$header{INFO}}, "${vt_string}" if (!$seen{vt});
	push @{$header{other}}, $command if ($command);

	foreach my $tag (@headorder) {
		foreach (@{$header{$tag}}) {
			next if (!$dpbounds && ($_ =~ /FILTER=<ID=LowDP/ || $_ =~ /FILTER=<ID=HighDP/));
			push @headerlines, $_;
		}
	}
	push @headerlines, $line;
	
	return @headerlines;
}
