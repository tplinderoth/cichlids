#!/usr/bin/perl

# extractVariants.pl

use warnings;
use strict;
use IO::Zlib;
use Getopt::Long;

my $version = '1.0.0';

my $skip_af;

die(qq/
extractVariants.pl v$version

extractVariants.pl [options] <input vcf>
OR
cat VCF | extractVariants.pl [options]

Input options:
--skip_af   Do not check that at least one allele has INFO\/AF in the interval \(0,1\) for retaining site

Notes:
* Assumes input VCF indels are left-aligned and normalized (bcftools norm -f <fasta>)
\n/) if (!@ARGV && -t STDIN);

GetOptions('--skip_af' => \$skip_af);
my %userargs = ('--skip_af' => $skip_af);

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

# declare variables to hold stats
my $nsites = 0; # number of sites in VCF
my $nsnps = 0; # total number SNPs
my $nsnp_sites; # total number sites containing SNPs
my $nsnp_widel = 0; # number of SNP sites within a deletion
my $nbisnps = 0; # number sites with biallelic snps
my $nbisnps_widel = 0; # number sites with biallelic snps within a deletion
my $ninsert = 0; # number of insertions
my $ndelete = 0; # number of deletions
my $nout_sites = 0; # number SNPs in output VCF

# proces VCF header

my $datestr = localtime();
my $command = "##extractVariantsCommand=<ID=extractVariants.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
        $command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

# print pretty header

my @header = makeHeader($vcffh, $command);

foreach my $hline (@header) {
	print STDOUT $hline;
}

# extract variants
my @indel = (0,0); # start and stop coordinates for indel
my ($end, $aa, $vartype, $chr);

$" = "\t";

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
        next if $tok[4] eq '.'; # skip nonvariable sites
        my @altarr = split(/,/, $tok[4]);
        my $issnp = 0;
	my $isins = 0;
	my $isdel = 0;
        my $altn;
        my @snp_alleles = ();
        my $refbase = substr $tok[3], 0, 1;
        my $i = 1;

        foreach my $a (@altarr) {
		if ($a eq '*') {
                        # do not record star alleles
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
			if ($altlen > $reflen) {
				$isins = 1;
				$ninsert++;
			} else {
				$isdel = 1;
				$ndelete++;
			}
                }
                $i++;
        }

        $nsnp_sites++ if ($issnp);
        $nsnp_widel++ if ($issnp && $widel);
	if (scalar(@snp_alleles) == 1) {
		$nbisnps++;
		$nbisnps_widel++ if ($widel);
	}

        # annotate INFO variant type
 	my $vt = '';
        if ($tok[7] =~ /VT=([^\s|;|=]+)/) {
         	$vt = $1;
		my $vt_original = $vt;
		$vt .= ',snp' if ($issnp && $vt !~ /snp/i);
                $vt .= ',deletion' if ($isdel && $vt !~ /deletion/i);
		$vt .= ',insertion' if ($isins && $vt !~ /insertion/i);
		$vt .= ',widel' if ($widel && $issnp && $vt !~ /widel/i);
		$tok[7] =~ s/VT=$vt_original/VT=$vt/;
	} else {
		$vt .= ',snp' if ($issnp);
                $vt .= ',deletion' if ($isdel);
                $vt .= ',insertion' if ($isins);
                $vt .= ',widel' if ($widel && $issnp);
		$vt =~ s/^,//;
		$tok[7] .= ";VT=$vt";
	}

	# print site
	my $isvar = 1;
	my $af = $1 if ($tok[7] =~ /[;\s]+AF=([^;\s]+)/);
	if (!$af && !$skip_af) {
		print STDERR "WARNING: No INFO/AF found for @tok[0..7]\n";
	} else {
		if (!$skip_af) {
			# check  whether at least one variant has frequency in the interval (0,1)
			$isvar = 0;
			foreach my $freq (split(',', $af)) {
				if ($freq > 0.0 && $freq < 1.0) {
					$isvar = 1;
					last;
				}
			}
		
		}
	}

	if ($isvar) {
		print STDOUT "@tok\n";
		$nout_sites++;
	}

}

close $vcffh;

my $nvariants = $nsnps + $ndelete + $ninsert; # total number of variants

# print stats

print STDERR "Last contig processed: $chr\n";
print STDERR "Number input VCF entries (lines without header): $nsites\n";
print STDERR "Number variants: $nvariants\n";
print STDERR "Number SNPs: $nsnps\n";
print STDERR "Number insertions: $ninsert\n";
print STDERR "Number deletions: $ndelete\n";
print STDERR "Number SNP sites: $nsnp_sites\n";
print STDERR "Number bialleleic SNP sites: $nbisnps\n";
print STDERR "Number SNP sites within a deletion: $nsnp_widel\n";
print STDERR "Number biallelic SNP sites within a deletion: $nbisnps_widel\n";
print STDERR "Number sites in output variant VCF: $nout_sites\n";

exit;

sub makeHeader {
        my ($fh, $command) = @_;

        # 1000 Genomes style VT annotation: ##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
        # widel denotes that the SNP is within a deletion
        my $vt_string = "##INFO=<ID=VT,Number=.,Type=String,Description=\"Denotes variant type\">\n";

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
                        push @headerlines, $_;
                }
        }
        push @headerlines, $line;
        
        return @headerlines;
}
