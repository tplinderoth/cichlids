#!/usr/bin/perl

# insertAnnotations.pl

use warnings;
use strict;
use Getopt::Long;
use IO::Zlib;

my $version = '0.0.2';
my $alfile = undef;
my $alfields = undef;
my $dpbounds = undef;
my $vcf = undef;
my $rmfilter = undef;

die(qq/
insertAnnotations.pl v$version

insertAnnotations.pl [options] <vcf file>
OR
cat <vcf file> | insertAnnotations.pl [options]

--alfile    File of outgroup alleles in pafAlleles format
--alfields  Column numbers of alleles to use from alleles file, 1-indexed and ','-separated
--dpbounds  LowDP and HighDP bounds to annotate FILTER filed, format '<INT lower bound>,<INT upper bound>'
--rmfilter  ','-delimited list of FILTER annotations to remove

Notes:
Assumes the same contig order in VCF and allele files
\n/) if (!@ARGV && -t STDIN);

GetOptions('alfile=s' => \$alfile, 'alfields=s' => \$alfields, 'dpbounds=s' => \$dpbounds, 'rmfilter=s' => \$rmfilter);

$vcf = pop @ARGV;

my %userargs = ('--alfile' => $alfile, '--alfields' => $alfields, '--dpbounds' => $dpbounds, '--rmfilter' => $rmfilter);

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

open(my $alfh, '<', $alfile) or die("Could not open allele file $alfile: $!\n");

my @alcol = map {$_ - 1} split(/,/, $alfields);
my $nalleles = scalar(@alcol);

my @dpcutoff = split(/,/, $dpbounds);
if (scalar(@dpcutoff) != 2) {
	die("ERROR: --dpbounds requires 2 values, a lower and an upper coverage threshold\n");
} elsif ($dpcutoff[0] >= $dpcutoff[1]) {
	die("ERROR: --dpbounds upper bound must be greater than the lower bound value\n");
}

my %delfilter;
if ($rmfilter) {
	foreach my $filter_anno (split(/,/,$rmfilter)) {
		$delfilter{$filter_anno} = 1;
	}
}

# print pretty VCF header

my $datestr = localtime();
my $command = "##insertAnnotationsCommand=<ID=insertAnnotations.pl,Version=${version},Date=\"$datestr\",CommandLineOptions=\"";
foreach my $arg (keys %userargs) {
	$command .= "$arg $userargs{$arg} " if (defined $userargs{$arg});
}
$command .= "$vcf" if ($vcf);
$command =~ s/\s$//;
$command .= "\">\n";

my $aa_string = "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n";
my $lowdp_string = "##FILTER=<ID=LowDP,Description=\"Site DP is less than median site DP -25%\">\n";
my $highdp_string = "##FILTER=<ID=HighDP,Description=\"Site DP is greater than median site DP +25%\">\n";

my @headorder = ('fileformat', 'reference', 'contig', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'other'); # header order

my %header = (fileformat => undef, ALT => undef, FILTER => undef, INFO => undef, FORMAT => undef, reference => undef, contig => undef, other => undef);

my %seen = (aa => 0, lowdp => 0, highdp => 0);

my $hline;
while (($hline = <$vcffh>) =~ /^##/) {
	if ($hline =~ /^##([^=]+)/i) {
		my $annotation = $1;
		
		if ($rmfilter && $annotation eq 'FILTER') {
			my $filter_id = $1 if ($hline =~ /^##FILTER=<ID=([^,]+)/);
			next if (exists $delfilter{$filter_id});
		}

		if ($alfile && $hline =~ /INFO=<ID=AA/) {
			# update AA annotation to ensure it is correct
			push @{$header{$annotation}}, $aa_string;
			++$seen{aa};
		} elsif ($dpbounds && $hline =~ /FILTER=<ID=LowDP/) {
			# update LowDP annotation
			push @{$header{$annotation}}, $lowdp_string;
			++$seen{lowdp};
		} elsif ($dpbounds && $hline =~ /FILTER=<ID=HighDP/) {
			# update HighDP annotation
			push @{$header{$annotation}}, $highdp_string;
			++$seen{highdp};
		} elsif ($annotation =~ /^fileformat|ALT|FILTER|INFO|FORMAT|reference|contig/) {
			push @{$header{$annotation}}, $hline;
		} else {
			push @{$header{other}}, $hline;
		}
	}
}

push @{$header{INFO}}, "${aa_string}" if ($alfile && !$seen{aa});
push @{$header{FILTER}}, "${lowdp_string}" if ($dpbounds && !$seen{lowdp});
push @{$header{FILTER}}, "${highdp_string}" if ($dpbounds && !$seen{highdp});
push @{$header{other}}, $command;

foreach my $tag (@headorder) {
	foreach (@{$header{$tag}}) {
		print STDOUT $_;
	}
}

print STDOUT $hline;

# process VCF sites

my $site_count = 0;
my $aa_count = 0;
my $low_count = 0;
my $high_count = 0;

$" = "\t";

<$alfh> if $alfile; # skip header
my $l1 = <$alfh>; # first site of new chr
die("No sites in allele file $alfile\n") if (!$l1);
my ($l2, $chr, $nsites) = (undef, "", 0);
my @outal;

while (<$vcffh>) {
	my @tok = split(/\s+/, $_);
	$site_count++;
	
	# ancestral allele
	if ($alfile) {
		if ($chr ne $tok[0]) {
		# store ancestral alleles in massive matrix
			@outal = ();
			$chr = $1 if ($l1 =~ /^(\S+)/);

			while ($chr ne $tok[0]) {
				my $l1 = <$alfh>;
				die("$tok[0] not found in allele file $alfile\n") if (!$l1);
				$chr = $1 if ($l1 =~ /^(\S+)/);
			}
			$l2 = $l1;
			my $c;

			do {
				push @outal, [(split(/\s+/, $l2))[@alcol]];
				if (defined($l2 = <$alfh>)) {
					$c = $1 if ($l2 =~ /^(\S+)/);
				}
			} while (defined($l2) && $c eq $chr);

			$l1 = $l2;
			$nsites = scalar(@outal);

			# debug
			# print "$nsites\n";
			#foreach (@outal) {
			#	print "@$_\n";
			#}
			#exit;
		}

		my %allele_counts = ();

		# collect reference alleles - expect only one
		die("More than one reference allele at $tok[0] $tok[1]: $tok[3]\n") if ($tok[3] =~ /,/);
		my $reflen = length($tok[3]);
		$allele_counts{$tok[3]} = 0;

		# collect alternate alleles
		foreach my $alt (split(/,/, $tok[4])) {
			$allele_counts{$alt} = 0;
		}

		# count number of times each ref and alt allele is represented among outgroups
		my $aa_seen = 0; # records whether at least one ref or alt allele is present among outgroups
		my $start = $tok[1]-1;
		my $end = $start + $reflen;

		# loop over each outgroup assembly
		for (my $assem = 0; $assem < $nalleles; $assem++) {
			
			# find counts for non-insertions
			my $standard_allele = ""; # collect alleles spanning ref length discarding deleted bases
			for (my $s = $start; $s < $end; $s++) {
				$standard_allele .= $outal[$s]->[$assem] if ($outal[$s]->[$assem] ne '*');
			}
			$standard_allele = '*' if (!$standard_allele);

			# check if standardized outgroup allele is among the ref and alt alleles
			if (exists $allele_counts{$standard_allele}) {
				$allele_counts{$standard_allele}++;
				$aa_seen = 1;
			}

			# check if outgroup allele is among alt insertions
			if ( length($outal[$start]->[$assem]) > $reflen) {
			# this is an outgroup insertion wrt to ref and so may be among alt alleles
				if (exists $allele_counts{$outal[$start]->[$assem]}) {
					$aa_seen = 1;
					$allele_counts{$outal[$start]->[$assem]}++;
				}
			}
		}

		my $anc_allele;
		if ($aa_seen) {
			my $maxcount = 0;
			my $nmax = 0;
			foreach my $allele (keys %allele_counts) {
				if ($allele_counts{$allele} >= $maxcount) {
					$nmax = $allele_counts{$allele} == $maxcount ? $nmax + 1:1;
					$maxcount = $allele_counts{$allele};
					$anc_allele = $allele;
				}
			}
			$anc_allele = "N" if ($nmax > 1 || $maxcount < 1);
		
		} else {
			$anc_allele="N";
		}

		# add ancestral allele to INFO
		if (!($tok[7] =~ s/AA=[^;|\s]+/AA=$anc_allele/)) {
			$tok[7] .= ";AA=$anc_allele";
			$aa_count++ if ($anc_allele ne 'N');
		}
	}

	# remove filter annotations
	if ($rmfilter) {
		if ($tok[6] ne '.') {
			my $filter_str = "";
			foreach my $fid (split(/;/, $tok[6])) {
				$filter_str .= "$fid;" unless exists $delfilter{$fid};
			}
			$filter_str =~ s/;*$//;
			$tok[6] = (!$filter_str) ? '.' : $filter_str; 
		}
	}

	# coverage filters
	if ($dpbounds) {
		my $dp = $1 if ($tok[7] =~ /DP=(\d+)/);
		if (defined $dp) {
			my $dptag = undef;
			if ($dp < $dpcutoff[0]) {
				$dptag = "LowDP";
				$low_count++;
			} elsif ($dp > $dpcutoff[1]) {
				$dptag = "HighDP";
				$high_count++;
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
		} else {
			print STDERR "No DP info for $tok[0] $tok[1], INFO=<$tok[7]>\n"
		}
	}

	print STDOUT "@tok\n";
}

close $vcffh;
close $alfh;

# print counts
print STDERR "Last contig processed: $chr\n";
print STDERR "Number sites in VCF: $site_count\n";
print STDERR "Number of sites with ancestral alleles: $aa_count\n" if ($alfile);

if ($dpbounds) {
	print STDERR "Number LowDP: $low_count\n";
	print STDERR "Number HighDP: $high_count\n";
}

exit;
