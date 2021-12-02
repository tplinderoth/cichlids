#!/usr/bin/perl

# variant_summary.pl

use warnings;
use strict;

die(qq/
variant_summary.pl <log directory> <log file prefix, '_[1-20,22-25].err' appended>
\n/) if (!@ARGV && scalar(@ARGV) < 2);

my $dir = $ARGV[0];
$dir =~ s/\/$//;
my $prefix = $ARGV[1];

my $nfields = 10;

my @counts;
# 0: Number of VCF entries (lines without header)
# 1: Number of variants (SNPs + indels)
# 2: Number insertions
# 3: Number deletions
# 4: Number SNPs
# 5: Number sites with at least one SNP
# 6: Number biallelic SNP sites
# 7: Number of deleted sites with at least one SNP
# 8: Number of biallelic SNP sites within a deletion
# 9: Number sites in biallelic VCF (biallelic SNPs outside deletions with INFO/AF in range (0,1))

print STDOUT qq/## CHR: Chromosome or scaffold
## ALL_SITES: Number of lines in all-sites VCF discounting header
## VARIANTS: Number of variants (SNPs + indels)
## INSERTIONS: Number of insertions
## DELETIONS: Number of deletions
## SNPS: Number single nucleotide changes
## SNP_SITES: Number sites with at least one SNP
## BISNP: Number biallelic SNP sites
## SNP_DEL: Number of sites overlapping a deletion with at least one SNP
## BISNP_DEL: Number of biallelic SNP sites within a deletion
## BIVCF: Number of sites in biallelic VCF (biallelic SNPs outside of deletions and with INFO\/AF in range (0,1))
\n/;

my @id = ('ALL_SITES', 'VARIANTS', 'INSERTIONS', 'DELETIONS', 'SNPS', 'SNP_SITES', 'BISNP', 'SNP_DEL', 'BISNP_DEL', 'BIVCF');

my @genome_total;
for (my $k = 0; $k < $nfields; $k++) {push @genome_total, 0;}

my %chr;

for (my $i=1; $i<26; $i++) {
	next if $i == 21;
	my @c;
	for (my $j=0; $j<$nfields; $j++) {push @c, 0;}
	my $logfile = "${dir}/${prefix}_${i}.err";
	if (!-f $logfile) {
		print STDERR "WARNING: $logfile does not exist\n";
		push @counts, \@c;
		next;
	} else {
		$chr{$i} = 1;
	}
	open(my $logfh, '<', $logfile) or die("Unable to open log file $logfile : $!\n");
	while (my $line = <$logfh>) {
		chomp $line;
		if ($line =~ /^Number VCF entries[^:]+:\s+(\d+)/) {
			$c[0] = $1;
		} elsif ($line =~ /^Number variants:\s+(\d+)/) {
			$c[1] = $1;
		} elsif ($line =~ /^Number insertions:\s+(\d+)/) {
			$c[2] = $1;
		} elsif ($line =~ /^Number deletions:\s+(\d+)/) {
			$c[3] = $1;
		} elsif ($line =~ /^Number SNPs:\s+(\d+)/) {
			$c[4] = $1;
		} elsif ($line =~ /^Number SNP sites:\s+(\d+)/) {
			$c[5] = $1;
		} elsif ($line =~ /^Number bialleleic SNP sites:\s+(\d+)/) {
			$c[6] = $1;
		} elsif ($line =~ /^Number SNP sites within a deletion:\s+(\d+)/) {
			$c[7] = $1;
		} elsif ($line =~ /^Number biallelic SNP sites within a deletion:\s+(\d+)/) {
			$c[8] = $1;
		} elsif ($line =~ /^Number sites in biallelic SNP VCF:\s+(\d+)/) {
			$c[9] = $1;
		}
	}
	close $logfh;
	if ($i < 24) {
		for (my $n = 0; $n < $nfields; $n++) {$genome_total[$n] += $c[$n];}
	}
	push @counts, \@c;
}

my $fmtstr = '%-17s';
for (my $n=0; $n<$nfields; $n++) {
	my $maxlen = length($id[$n]);
	my @len = (length($counts[22][$n]), length($counts[23][$n]), length($genome_total[$n]));
	foreach my $l (@len) { $maxlen = $l if $l > $maxlen;}
	$maxlen += 3;
	$fmtstr .= "%-${maxlen}s"
}

printf("${fmtstr}\n", "CHR", @id);
my $chrn = 0;
for (my $i=0; $i<22; $i++) {
	$chrn += $chrn != 20 ? 1:2;
	printf("${fmtstr}\n", "chr${chrn}", @{$counts[$i]}) if exists $chr{$chrn};
}
printf("${fmtstr}\n", "Genome_main_chr", @genome_total);
printf("${fmtstr}\n", "U_scaffolds", @{$counts[22]}) if exists $chr{24};
printf("${fmtstr}\n", "chrM", @{$counts[23]}) if exists $chr{25};

exit;
