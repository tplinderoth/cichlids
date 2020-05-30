#!/usr/bin/perl

# queryMeta.pl

# requires bcftools (installed and in PATH)

use warnings;
use strict;

my $version = '1.0.0';

die(qq/
callsetTools.pl $version 

callsetTools.pl [command]

commands:
sampleID     retrieve Sanger IDs for samples belonging to a group
alleleFreq   calculate allele frequencies 
\n/) if (!@ARGV || $ARGV[0] eq '--help' || $ARGV[0] eq '-h');


my $command = shift @ARGV;

if ($command eq 'sampleID') {
	sampleid();
} elsif ($command eq 'alleleFreq') {
	alleleFreq();
} else {
	die("Unknown command $command\n");
}

exit;

sub alleleFreq {

die(qq/
callsetTools.pl alleleFreq [chr number] [position] [sample list] [vcf]

Sample list:
2-column file in sampleID output format where each row lists (1) Sanger ID of sample, (2) group identity.
Assumes a header line is present. Allele frequencies will be calculated for each group.

VCF:
Indexed VCF file.
\n/) if (!@ARGV || scalar(@ARGV) < 4);

my ($chr, $pos, $idlist, $vcf) = @ARGV;

die("VCF file $vcf does not exist\n") if (!-e $vcf);

my %groups;
my %seen;
my @order;
open(my $listfh, '<', $idlist) or die("Couldn't open sample list $idlist: $!\n");
<$listfh>; #skip header
while (<$listfh>) {
	chomp;
	if ($_ =~ /^(\S+)\s+(\S+)/) {
		$groups{$1} = $2;
		if (!exists $seen{$2}) {
			push @order, $2;
			$seen{$2} = 1;
		}
	}
}

chomp(my $head = `bcftools view -h $vcf | tail -n -1`);
my %groupidx;
my $i = 0;
foreach (split(/\s+/, $head)) {
	push @{$groupidx{$groups{$_}}}, $i if (exists $groups{$_});
	$i++;
}

foreach(keys %groupidx) {
	@{$groupidx{$_}} = sort {$a <=> $b} @{$groupidx{$_}};
	#print "$_: @{$groupidx{$_}}\n"; # debug
}

my @vcfpos = split(/\s+/,`bcftools view -H -r chr$chr:$pos $vcf | tail -n 1`);

my $aa = 'N';
$aa = $1 if ($vcfpos[7] =~ /AA=([^;|\s]+)/);
print STDOUT "$vcfpos[0]\t$vcfpos[1]\tref=$vcfpos[3]\talt=$vcfpos[4]\tanc=$aa\t$vcfpos[6]\t$vcfpos[7]\n";

print STDOUT "group\talt_freq\tn\tndata\n";

foreach my $id (@order) {
	#print "$id: @{$groupidx{$id}}\n"; # debug
	my $n = 0;
	my $ndata = 0;
	my $altsum = 0;
	foreach my $ind (@vcfpos[@{$groupidx{$id}}]) {
		if ($ind =~ /^([0|1]{1}[\/|\|]{1}[0|1]{1})/) {
			my $geno = $1;
			my $nalt = () = $geno =~ /1/g;
			$altsum += $nalt;
			$ndata++;
		}
		$n++;
	}
	my $altfreq = $altsum / (2*$ndata);
	print STDOUT "$id\t$altfreq\t$n\t$ndata\n";
}

}

sub sampleid {

die(qq/
callsetTools.pl sampleID [type] [list] [metadata file]

Types:
species   return IDs of samples belonging to a species specified like 'Genus Species' (must include single-quotes)

List:
Either a space-delimited list of names or 1-column file listing 1 name per row.
Names with spaces must be enclosed in single quotes if not using a file list, e.g. 'Genus\\s+Species'.
Quotes from names that use them can be excluded, e.g. 
'Maylandia 'lanisticola north'' and 'Maylandia lanisticola north' are both accepted.

Metadata file:
Tab-delimited file of metadata.
\n/) if (!@ARGV || scalar (@ARGV) < 3);

my $type = shift @ARGV;

my $metafile = pop @ARGV;
open(my $metafh, '<', $metafile) or die("Couldn't open metadata file $metafile: $!\n");

my $fname;
my %names;
my @order;

if (-f $ARGV[0]) {
	open(my $fh, '<', $ARGV[0]) or die("Couldn't open names file $ARGV[0]: $!\n");
	while (<$fh>) {
		chomp;
		$_ =~ s/\s+/_/g;
		$_ =~ s/'|"//g;
		push @order, $_;
		@{$names{$_}} = ();
	}
	close $fh;
} else {
	foreach (@ARGV) {
		$_ =~ s/\s+/_/g;
                $_ =~ s/'|"//g;
		push @order, $_;
		@{$names{$_}} = ();
	}
}

#foreach (@order) {print "$_\n";}; exit; # debug

my %metahead;
my $c = 0;
foreach my $head (split(/\s+/, <$metafh>)) {
	$metahead{$head} = $c;
	$c++;
}

my @idx;
if ($type eq 'species') {
	@idx = ($metahead{genus}, $metahead{species});
	print STDOUT "primary_id\tspecies\n";
} else {
	die("Invalid type $type\n");
}

while (<$metafh>) {
	my @l = split("\t", $_);
	my $name = join('_', @l[@idx]);
	$name =~ s/\s+/_/g;
	$name =~ s/'|"//g;
	push @{$names{$name}}, $l[$metahead{primary_id}];
}
close $metafh;

foreach my $groupid (@order) {
	foreach my $sampid (@{$names{$groupid}}) {
		print STDOUT "$sampid\t$groupid\n";
	}
}

}
