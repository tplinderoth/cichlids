#!/usr/bin/perl

# annotate_summary.pl

use warnings;
use strict;

die(qq/
annotate_summary.pl <log directory> <log file prefix, '_[1-20,22-25].err' appended>
\n/) if (!@ARGV && scalar(@ARGV) < 2);

my $dir = $ARGV[0];
$dir =~ s/\/$//;
my $prefix = $ARGV[1];

my ($lines_total, $fail_total, $aa_total) = (0,0,0);
my ($lines_total_u, $fail_total_u, $aa_total_u) = (0,0,0);
my ($lines_total_mt, $fail_total_mt, $aa_total_mt) = (0,0,0);
my ($aaperc, $npass, $passperc, $id) = (0, 0, 0, '');
my (%filters, %filters_u, %filters_mt);
my (%flags, %flags_u, %flags_mt);
my ($flags_total, $flags_total_u, $flags_total_mt) = (0,0,0);
my @l;

for (my $i=1; $i < 26; $i++) {
	next if $i == 21;
	my $logfile = "${dir}/${prefix}_${i}.err";
	if (!-f $logfile) {
		print STDERR "WARNING: $logfile does not exist\n";
		next;
	}
	open(my $logfh, '<', $logfile) or die("Unable to open log file $logfile : $!\n");
	my ($nlines, $naa, $nfail, $configs, $flags) = (0,0,0,0,0);
	while (<$logfh>) {
		chomp;
		my $line = $_;
		$nlines = $1 if ($line =~ /^Number sites in VCF:\s+(\d+)/);
		$naa = $1 if ($line =~ /^Number of sites with ancestral alleles:\s(\d+)/);
		$nfail = $1 if ($line =~ /Number of FAIL sites:\s(\d+)/);
		if ($line =~ /^===== FILTER CONFIGURATION COUNTS =====/) {
			$configs = 1;
			next;
		}
		if ($line =~ /^===== FILTER FLAG COUNTS =====/) {
			$flags = 1;
			next;
		}
		if ($configs == 1 && $flags == 0 && $line =~ /\d+$/) {
			@l = split(/\s+/, $line);
			$l[0] =~ s/:$//;
			my @f = split(';', $l[0]);
			foreach my $fid (@f) {
				if ($i < 24) {
					if (exists $flags{$fid}) {
						$flags{$fid} += $l[1];
					} else {
						$flags{$fid} = $l[1];
					}
					$flags_total += $l[1];
				} elsif ($i == 24) {
					if (exists $flags_u{$fid}) {
                                        	$flags_u{$fid} += $l[1];
                                       	} else {
                                        	$flags_u{$fid} = $l[1];
                                       	}
					$flags_total_u += $l[1];
				} elsif ($i == 25) {
					 if (exists $flags_mt{$fid}) {
                                                $flags_mt{$fid} += $l[1];
                                       	} else {
                                                $flags_mt{$fid} = $l[1];
                                       	}
					$flags_total_mt += $l[1];
				} else {
					die("Unknown chromosome value, $i\n");
				}
			}
			my $fcomb = join(';',sort {$a cmp $b} @f);
			if ($i < 24) {
				$filters{$fcomb} = exists $filters{$fcomb} ? $filters{$fcomb} + $l[1] : $l[1];
			} elsif ($i == 24) {
				$filters_u{$fcomb} = exists $filters_u{$fcomb} ? $filters_u{$fcomb} + $l[1] : $l[1];
			} elsif ($i == 25) {
				$filters_mt{$fcomb} = exists $filters_mt{$fcomb} ? $filters_mt{$fcomb} + $l[1] : $l[1];
			} else {
				die("Unknown chromosome value, $i\n");
			}
		}
	}
	close $logfh;
	my $id = '';
	if ($i < 24) {
		$id = "chr${i}";
		$lines_total += $nlines;
		$aa_total += $naa;
		$fail_total += $nfail;
	} elsif ($i == 24) {
		$id = "U_scaffolds";
		$lines_total_u += $nlines;
                $aa_total_u += $naa;
                $fail_total_u += $nfail;
	} elsif ($i == 25) {
		$id = "chrM";
		$lines_total_mt += $nlines;
                $aa_total_mt += $naa;
                $fail_total_mt += $nfail;
	} else {
		die("Unknown chromosome value, $i\n");
	}
	$aaperc = $nlines > 0 ? sprintf("%.1f", ($naa/$nlines)*100) : 0;
	$npass = $nlines - $nfail;
	$passperc = $nlines > 0 ? sprintf("%.1f", ($npass/$nlines)*100) : 0;
	print STDOUT "$id:NLINES=$nlines;NAA=$naa(${aaperc}%);PASS=$npass(${passperc}%)\n" if ($i < 24);
}

if ($lines_total > 0) {
	$aaperc = sprintf("%.1f", ($aa_total/$lines_total)*100);
        $npass = $lines_total - $fail_total;
	$passperc = sprintf("%.1f", ($npass/$lines_total)*100);
	print STDOUT "main_chromosomes:NLINES=$lines_total;NAA=$aa_total(${aaperc}%);PASS=$npass(${passperc}%)\n";
}

if ($lines_total_u > 0) {
        $aaperc = sprintf("%.1f", ($aa_total_u/$lines_total_u)*100);
        $npass = $lines_total_u - $fail_total_u;
        $passperc = sprintf("%.1f", ($npass/$lines_total_u)*100);
	print STDOUT "U_scaffolds:NLINES=$lines_total_u;NAA=$aa_total_u(${aaperc}%);PASS=$npass(${passperc}%)\n";
}

if ($lines_total_mt > 0) {
        $aaperc = sprintf("%.1f", ($aa_total_mt/$lines_total_mt)*100);
        $npass = $lines_total_mt - $fail_total_mt;
        $passperc = sprintf("%.1f", ($npass/$lines_total_mt)*100);
        print STDOUT "chrM:NLINES=$lines_total_mt;NAA=$aa_total_mt(${aaperc}%);PASS=$npass(${passperc}%)\n";
}

print STDOUT "\n===== FILTER CONFIGURATION SUMMARY =====\n";

my $fperc;
if ($lines_total > 0) {
	print STDOUT "\n= MAIN CHROMOSOMES =\n";
	foreach (keys %filters) {
		$fperc = sprintf("%.1f", ($filters{$_}/$fail_total)*100);
		print STDOUT "$_: $filters{$_} (${fperc}%)\n";
	}
}

if ($lines_total_u > 0) {
        print STDOUT "\n= UNPLACED SCAFFOLDS =\n";
        foreach (keys %filters_u) {
                $fperc = sprintf("%.1f", ($filters_u{$_}/$fail_total_u)*100);
                print STDOUT "$_: $filters_u{$_} (${fperc}%)\n";
        }
}

if ($lines_total_mt > 0) {
        print STDOUT "\= MITOCHONDRIAL GENOME =\n";
        foreach (keys %filters_mt) {
                $fperc = sprintf("%.1f", ($filters_mt{$_}/$fail_total_mt)*100);
                print STDOUT "$_: $filters_mt{$_} (${fperc}%)\n";
        }
}

print STDOUT "\n===== FLAG COUNTS (among filter configurations) =====\n";

if ($lines_total > 0) {
        print STDOUT "\n= MAIN CHROMOSOMES =\n";
        foreach (keys %flags) {
                $fperc = sprintf("%.1f", ($flags{$_}/$flags_total)*100);
                print STDOUT "$_: $flags{$_} (${fperc}%)\n";
        }
}

if ($lines_total_u > 0) {
        print STDOUT "\n= UNPLACED SCAFFOLDS =\n";
        foreach (keys %flags_u) {
                $fperc = sprintf("%.1f", ($flags_u{$_}/$flags_total_u)*100);
                print STDOUT "$_: $flags_u{$_} (${fperc}%)\n";
        }
}

if ($lines_total_mt > 0) {
        print STDOUT "\n= MITOCHONDRIAL GENOME =\n";
        foreach (keys %flags_mt) {
                $fperc = sprintf("%.1f", ($flags_mt{$_}/$flags_total_mt)*100);
                print STDOUT "$_: $flags_mt{$_} (${fperc}%)\n";
        }
}


exit;
