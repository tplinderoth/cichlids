#!/usr/bin/perl

# insertAnnotations.pl

use warnings;
use strict;
use Getopt::Long;
use IO::Zlib;
#use warnings FATAL => 'all';

my $version = '0.2.0';
my $alfile = undef;
my $alfields = undef;
my $anctype = 'parsimony';
my $dpbounds = undef;
my $vcf = undef;
my $rmfilter = undef;
my $fstart = 4;
my $bedgraph = undef;
my $bedlist = undef;

die(qq/
insertAnnotations.pl version $version

insertAnnotations.pl [options] <vcf file>
OR
cat <vcf file> | insertAnnotations.pl [options]

--alfile    File of outgroup alleles in pafAlleles format
--alfields  Column numbers of alleles to use from alleles file, 1-indexed and ','-separated [$fstart...]
--anctype   Infer ancestral allele using 'parsimony' or as the 'major' allele among outgroups [$anctype]
--dpbounds  LowDP and HighDP bounds to annotate FILTER field, format '<INT lower bound>,<INT upper bound>'
--rmfilter  ','-delimited list of FILTER annotations to remove
--bedgraph  Bed format file with annotations to add to the FILTER field
--bedlist   ','-delimited list of FILTER annotations to use from bedgraph (ignore all others)

Notes:
*Assumes the same contig order in VCF and allele files.

*For inferring ancestral alleles with parsimony it is assumed that outgroups in the pafAlleles file 
 are ordered with increasing divergence from the ingroup (greater column number = greater divergence).
 Alternatively, the outgroup columns can be passed in order of their increasing divergence from the
 ingroup using the --alfields option.

*Parsimony currently only implemented for 3 outgroup case.

*Particular FILTER annotations will be updated if the function that applies them (--dpbounds, --bedgraph)
 is called.

*Bedgraph must be sorted in same order as VCF.
\n/) if (-t STDIN && !@ARGV);

GetOptions('alfile=s' => \$alfile, 'alfields=s' => \$alfields, 'dpbounds=s' => \$dpbounds, 'rmfilter=s' => \$rmfilter, 'anctype=s' => \$anctype,
'bedgraph=s' => \$bedgraph, 'bedlist=s' => \$bedlist);

# this is for printing the command to the VCF header
my %userargs = ('--alfile' => $alfile, '--alfields' => $alfields, '--dpbounds' => $dpbounds, '--rmfilter' => $rmfilter, '--anctype' => $anctype,
'--bedgraph' => $bedgraph, '--bedlist' => $bedlist);
# disable printing commands if not applicable
$userargs{"--anctype"} = undef if (!$alfile);

# read in VCF

$vcf = pop @ARGV;

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

# initalize ancestral allele inference

my ($alfh, $nalleles, $aamethod, @outal, @col, @parstable, %tree_counts);
my ($l1, $l2, $chr, $nsites) = (undef, undef, "", 0);

if ($alfile) {
	open($alfh, '<', $alfile) or die("Could not open allele file $alfile: $!\n");

	if ($alfields) {
		@col = map {die("Argument to --alfields must be ','-delimited integers > 0") if $_ =~ /\D/; $_ - 1} split(/,/, $alfields);
                $nalleles = scalar @col;
	} else {
		$nalleles = split(/\s+/, <$alfh>)-3;
		@col = ($fstart-1) .. ($fstart + $nalleles - 2);
	}

        $l1 = <$alfh>; # first site of new chr
        die("No sites in allele file $alfile\n") if (!$l1);

	if (lc($anctype) eq 'parsimony') {
		$aamethod = 1;
		@parstable = makeLookup(\%tree_counts);
		treeCounter(\%tree_counts, $nalleles, 0);
	} elsif (lc($anctype) eq 'major') {
		$aamethod = 2;
	} else {
		die("Invalid --anctype argument: $anctype\n");
	}
} else {
	$anctype = undef;
}

# initialize filter annotation

my (@dpcutoff, %delfilter, $bedfh);

if ($dpbounds) {
	
	@dpcutoff = split(/,/, $dpbounds);
	if (scalar(@dpcutoff) != 2) {
		die("ERROR: --dpbounds requires 2 values, a lower and an upper coverage threshold\n");
	} elsif ($dpcutoff[0] >= $dpcutoff[1]) {
		die("ERROR: --dpbounds upper bound must be greater than the lower bound value\n");
	}
	
	$delfilter{HighDP} = 1;
	$delfilter{LowDP} = 1;
}

my %bedfilter;
if ($bedgraph) {
	die("Must provide --bedlist if using a bedgraph file of filter annotations\n") if (!$bedlist);
	open($bedfh, '<', $bedgraph) or die("Unable to open bedgraph file $bedgraph: $!\n");
	foreach (split(/,/, $bedlist)) {
		$delfilter{$_} = 1;
		$bedfilter{$_} = 1;
	}
	$bedlist =~ s/,/\|/g;
}

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
my $lowcov_string = "##FILTER=<ID=LowCov,Description=\"Total coverage among QC subset individuals below 0.05 quantile\">\n";
my $excesscov_string = "##FILTER=<ID=ExcessCov,Description=\"Total coverage among QC subset individuals above 0.95 quantile\">\n";
my $lowmq_string = "##FILTER=<ID=LowMQ,Description=\"Root mean square map quality below 35 for QC subset individuals\">\n";
my $mqzero_string = "##FILTER=<ID=ExcessMQ0,Description=\"Fraction of reads with map quality zero among QC subset individuals greater than 0.2\">\n";
my $nocov_string = "##FILTER=<ID=NoCov,Description=\"No mapped reads among QC subset individuals\">\n";

my @headorder = ('fileformat', 'reference', 'contig', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'other'); # header order

my %header = (fileformat => undef, ALT => undef, FILTER => undef, INFO => undef, FORMAT => undef, reference => undef, contig => undef, other => undef);

my %seen = (aa => 0, lowdp => 0, highdp => 0);

my $hline;
my $filter;
while (($hline = <$vcffh>) =~ /^##/) {
	if ($hline =~ /^##([^=]+)/i) {
		my $annotation = $1;

		if (($rmfilter || $bedgraph) && $hline =~ /^##FILTER=<ID=([^,]+)/) {
			next if exists $delfilter{$1}
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

push @{$header{INFO}}, $aa_string if ($alfile && !$seen{aa});
push @{$header{FILTER}}, $lowdp_string if ($dpbounds && !$seen{lowdp});
push @{$header{FILTER}}, $highdp_string if ($dpbounds && !$seen{highdp});
push @{$header{FILTER}}, $nocov_string if ($bedfilter{NoCov});
push @{$header{FILTER}}, $lowcov_string if ($bedfilter{LowCov});
push @{$header{FILTER}}, $excesscov_string if ($bedfilter{ExcessCov});
push @{$header{FILTER}}, $lowmq_string if ($bedfilter{LowMQ});
push @{$header{FILTER}}, $mqzero_string if ($bedfilter{ExcessMQ0});

push @{$header{other}}, $command;

foreach my $tag (@headorder) {
	foreach (@{$header{$tag}}) {
		print STDOUT $_;
	}
}

print STDOUT $hline;

# process VCF sites

my %filtercounts;
map { $filtercounts{$1} = 0 if $_ =~ /ID=([^,]+)/ } @{$header{FILTER}};

my $vcfchr;
my @bedtok = (0,0,0);
my $site_count = 0;
my $aa_count = 0;
my @tok;

$" = "\t";

while (<$vcffh>) {
	@tok = split(/\s+/, $_);
	$site_count++;

	# update begraph region
	if ($bedgraph) {
		$vcfchr = $1 if ($tok[0] =~ /(\d+)$/);
		while (!eof($bedfh) && $bedtok[0] < $vcfchr) {
			bedsplit(\@bedtok, $bedfh);
		}
		while (!eof($bedfh) && $bedtok[2] < $tok[1]) {
			bedsplit(\@bedtok, $bedfh);
			last if ($bedtok[0] != $vcfchr);
		}
	}
	
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
				push @outal, [(split(/\s+/, $l2))[@col]];
				if (defined($l2 = <$alfh>)) {
					$c = $1 if ($l2 =~ /^(\S+)/);
				}
			} while (defined($l2) && $c eq $chr);

			$l1 = $l2;
			$nsites = scalar(@outal);

			# debug
			#print "$nsites\n";
			#foreach (@outal) {
			#	print "@$_\n";
			#}
			#exit;
		}

		# collect ingroup alleles
		die("More than one reference allele at $tok[0] $tok[1]: $tok[3]\n") if ($tok[3] =~ /,/);
		my $reflen = length($tok[3]);
		my @ingroup = ($tok[3], split(/,/, $tok[4]));

		# collect alleles for each outgroup over the length of the ingroup ref allele
		my @outgroup;
		my $start = $tok[1]-1;
		my $end = $start + $reflen;
		# loop over each outgroup assembly
		for (my $assem=0; $assem<$nalleles; $assem++) {

			if (length($outal[$start]->[$assem]) > $reflen) {
				# insert wrt to ingroup reference
				push @outgroup, $outal[$start]->[$assem];
			} else {
				# non-insertion wrt to ingroup reference
				$outgroup[$assem] = "";
				for (my $s=$start; $s<$end; $s++) {
					$outgroup[$assem] .= $outal[$s]->[$assem] if ($outal[$s]->[$assem] ne '*');
				}

				if (!$outgroup[$assem]) {
					# check if all aligned sites were deleted
					$outgroup[$assem] = '*';
				} else {
					# check for deletions longer that ingroup deletions
					if ($outal[$end-1]->[$assem] eq '*' && $end <= $#outal && $outal[$end]->[$assem] eq '*') {
						$outgroup[$assem] .= '*'; # denotes that next site beyond last reference allele site is deleted
					}
				}
			}
		}

		# find ancestral allele
		my $anc_allele;

		if ($aamethod == 1) {
			# use parsimony to find ancestral allele

			# Find the outgroup matching configuration for each ingroup allele in the lookup table
			my @outsp;
			my $tree;

			foreach my $inallele (@ingroup) {

				my $tabref = \@parstable;
				$tree = "";

				foreach my $outallele (@outgroup) {
					if ($outallele =~ /N/ ) {
						# outgroup allele is missing
						$tabref = ${$tabref}[2];
						$tree .= 2;
					} else {
						# check if ingroup allele matches outgroup
						my $match = $inallele eq $outallele ? 1 : 0;
						$tabref = ${$tabref}[$match];
						$tree .= $match;
					}

				}

				push @outsp, $tabref;
				$tree_counts{$tree}++;
				$tree .= ',';
			}

			# resolve exceptional cases for most parsimonious solution
			if ($nalleles == 3) {
				# 3 outgroup case
				if ($tree =~ /100/ && $tree !~ /011/) {
					push @outsp, 0;
				}
			}

			# pick ancestral allele from the most related outgroup species
			my $spid = (sort {$a <=> $b} @outsp)[0];
			$anc_allele = $spid != 9 ? $outgroup[$spid] : 'N';

			# debug
			#print "$tok[0]\t$tok[1]\t@ingroup;\t@outgroup;\t@outsp;\t$anc_allele\n"; next;			

		} elsif ($aamethod == 2) {
			# use major allele among outgroups as ancestral allele

			my %allele_counts = ();
			foreach my $inallele (@ingroup) {
				$allele_counts{$inallele} = 0;
			}

			# count occurances of outgroup alleles among the ingroup
			my $aa_seen = 0; # records whether at least one ref or alt allele is present among outgroups
			foreach my $outallele (@outgroup) {
				if (exists $allele_counts{$outallele}) {
					$allele_counts{$outallele}++;
					$aa_seen = 1;				
				}
			}

			if ($aa_seen) {
				my $maxcount = 0;
				my $nmax = 0;
					foreach my $allele (@ingroup) {
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
		}

		# add ancestral allele to INFO
		if (!($tok[7] =~ s/AA=[^;|\s]+/AA=$anc_allele/)) {
			$tok[7] .= ";AA=$anc_allele";
		}
		$aa_count++ if ($anc_allele ne 'N');
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

	# apply filters
	if ($dpbounds || $bedgraph) {
		my @filter_annotations;
	
		# coverage filters
		if ($dpbounds) {
			my $dp = $1 if ($tok[7] =~ /DP=(\d+)/);
			if (defined $dp) {
				if ($dp < $dpcutoff[0]) {
					push @filter_annotations, "LowDP";
				} elsif ($dp > $dpcutoff[1]) {
					push @filter_annotations, "HighDP";
				}
			} else {
				print STDERR "No DP info for $tok[0] $tok[1], INFO=<$tok[7]>\n";
			}
		}

		# bedgraph filters
		if ($bedgraph && $vcfchr == $bedtok[0] && $tok[1] > $bedtok[1] && $tok[1] <= $bedtok[2]) {
			for (my $i=3; $i <= $#bedtok; $i++) {
				push @filter_annotations, $bedtok[$i] if ($bedtok[$i] ne "0" && exists $bedfilter{$bedtok[$i]});
			}
		}

		# update FILTER field with new annotations
		if (@filter_annotations) {
			if ($tok[6] eq "." || $tok[6] eq "PASS") {
				$tok[6] = join(';',@filter_annotations);
			} else {
				$tok[6] .= ";" . join(';',@filter_annotations);
			}
		} else {
			$tok[6] = "PASS" if ($tok[6] eq '.');
		}
	}
	map {$filtercounts{$_}++ unless $_ eq '.'} split(':', $tok[6]);

	print STDOUT "@tok\n";
}

close $vcffh;
close $alfh if ($alfile);
close $bedfh if ($bedgraph);

# print counts
print STDERR "Last contig processed: $tok[0]\n";
print STDERR "Number sites in VCF: $site_count\n";

if ($alfile) {
	print STDERR "Number of sites with ancestral alleles: $aa_count\n";
	if ($aamethod == 1) {
		print STDERR "Ingroup allele matching configurations\n";
		treeCounter(\%tree_counts, $nalleles, 1);
	}
}

map { print STDERR "Number $1: $filtercounts{$1}\n" if $_ =~ /ID=([^,]+)/ } @{$header{FILTER}};

exit;

sub bedsplit {
	my ($tok, $fh) = @_;
	@$tok = split(/\s+/, readline($$fh));
	$$tok[0] = $1 if $$tok[0] =~ /(\d+)$/;
}

sub makeLookup {
### returns species representing the ancestral allele ###

### species in descending relatedness to Astatotilapia calliptera
# Simochromis diagramma, Tropheini tribe (0)
# Neolamprologus brichardi, Lamprologini tribe (1)
# Oreochromis niloticus (2)
# not able to assign (-9)

### codes corresponding to elements of of the lookup table
# [S. diagramma]->[N. brichardi]->[O. niloticus]
# 0: mismatch to A. calliptera allele
# 1: match A. calliptera allele
# 9: not able to assign

        my @anc;

        $anc[0]->[0]->[0] = 9;
        $anc[0]->[0]->[1] = 9;
        $anc[0]->[0]->[2] = 9;
        $anc[0]->[1]->[0] = 9;
        $anc[0]->[1]->[1] = 1;
        $anc[0]->[1]->[2] = 1;
        $anc[0]->[2]->[0] = 9;
        $anc[0]->[2]->[1] = 2;
        $anc[0]->[2]->[2] = 9;
        $anc[1]->[0]->[0] = 9;
        $anc[1]->[0]->[1] = 0;
        $anc[1]->[0]->[2] = 0;
        $anc[1]->[1]->[0] = 0;
        $anc[1]->[1]->[1] = 0;
        $anc[1]->[1]->[2] = 0;
        $anc[1]->[2]->[0] = 0;
        $anc[1]->[2]->[1] = 0;
        $anc[1]->[2]->[2] = 0;
        $anc[2]->[0]->[0] = 9;
        $anc[2]->[0]->[1] = 2;
        $anc[2]->[0]->[2] = 9;
        $anc[2]->[1]->[0] = 1;
        $anc[2]->[1]->[1] = 1;
        $anc[2]->[1]->[2] = 1;
        $anc[2]->[2]->[0] = 9;
        $anc[2]->[2]->[1] = 2;
        $anc[2]->[2]->[2] = 9;

        if (@_) {
                my $config_counts = $_[0];
                for (my $i=0; $i<3; $i++) {
                        for (my $j=0; $j<3; $j++) {
                                for (my $k=0; $k<3; $k++) {
                                        $$config_counts{"$i$j$k"} = 0;
                                }
                        }
                }
        }

	return @anc;
}

sub treeCounter {
# fun = 0: initalize count hash
# fun = 1: print count hash

	my ($counts, $nout, $fun) = @_;

	# generate array of configurations
	my @config = ("");

	for (my $i=0; $i<$nout; $i++) {
		my @arr;
		foreach my $str (@config) {
			foreach my $state ((0,1,2)) {
				push @arr, "$str$state";
			}
		}
		@config = @arr;
	}

	# process configurations
	foreach (@config) {
		if ($fun == 0) {
			$$counts{$_} = 0;
		} elsif ($fun == 1) {
			print STDERR "$_: $$counts{$_}\n";
		}
	}
}
