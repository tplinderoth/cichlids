#!/usr/bin/perl

# queryMeta.pl

# requires bcftools and vcftools (installed and in PATH) for sum functions

use warnings;
use strict;
use Getopt::Long;
use Getopt::Std;

my $version = '1.2.0';

die(qq/
callsetTools.pl $version

Requires bcftools and vcftools to be installed and in \$PATH

callsetTools.pl [command]

commands:
sampleID     retrieve Sanger IDs for samples belonging to a group
alleleFreq   calculate allele frequencies
multiFreq    calculate allele frequencies for a list of sites
convert      convert VCF to other formats (allows subsetting)
\n/) if (!@ARGV || $ARGV[0] eq '--help' || $ARGV[0] eq '-h');


my $command = shift @ARGV;

if ($command eq 'sampleID') {
	sampleid();
} elsif ($command eq 'alleleFreq') {
	alleleFreq();
} elsif ($command eq 'multiFreq') {
	multiFreq();
} elsif ($command eq 'convert') {
	convert();
} else {
	die("Unknown command $command\n");
}

exit;

sub convert {

my ($subset, $pass, $mafcount, $scale, $chr) = (undef, undef, 0, undef, undef);

die(qq/
callsetTools convert [options] [format] [vcf]

Format:
beaglepl    Beagle genotype likelihood format
beagleDist  Modified Beagle genotype likelihood format for ngsDist
geno        ANGSD geno format

Options:
subset      File listing individuals to include
pass        Include only sites that 'PASS' according to FILTER field
mafcount    Discard sites with fewer than INT minor allele counts [$mafcount]
scale       Scale likelihoods such that they sum to 1 for each individual
chr         Comma-separated list of chromosome numbers (required for beaglepl and geno formats)

Notes:
*Keeps only biallelic sites.
\n/) if (!@ARGV || scalar(@ARGV) < 2);

my $vcf = pop @ARGV;
die ("Unable to locate VCF $vcf: $!\n") if (!-f $vcf);

my $format = pop @ARGV;

GetOptions('subset=s' => \$subset, 'pass' => \$pass, 'mafcount=i' => \$mafcount, 'scale' => \$scale, 'chr=s' => \$chr);
my @chrn = sort {$a <=> $b} split(',',$chr) if ($chr);

my $command="bcftools view -m2 -M2";
$command .= " -f PASS" if ($pass);
if ($subset) {
	die("Unable to locate file of subset IDs $subset\n") if (!-f $subset);
	$command .= " -S $subset -a";
}
$command .= " -c ${mafcount}:minor";
if ($chr) {
	$command .= " -r ";
	foreach my $num (@chrn) {
		$command .= "chr${num},"
	}
	$command =~ s/,$//;
}
$command .= " $vcf";

my $fmt;
if ($format eq 'beaglepl' || $format eq 'geno' || $format eq 'beagleDist') {
	$command .= " | vcftools --vcf - --BEAGLE-PL --stdout";
	if ($chr) {
		foreach my $num (@chrn) {
			$command .= " --chr chr$num";
		}
	} else {
		die("$format requires --chr\n");
	}
	$command .= " |";

	if ($format eq 'beaglepl') {
		$fmt = 1;
	} elsif ($format eq 'geno') {
		$fmt = 2;
	} elsif ($format eq 'beagleDist') {
		$fmt = 3;
	}

} else {
	die("Unrecognized format $format\n");
}

print STDERR "Reading from '$command'\n";

open(my $stream, $command) or die("Unable to read input stream: $!\n");

$" = "\t";
while (<$stream>) {
	if ($fmt == 1 ) {
		#BEAGLEPL
		print STDOUT $_;
	} elsif ($fmt == 2 || $fmt == 3) {
		# GENO OR NGSDIST BEAGLEPL
		next if ($_ =~ /^marker\b/);
		my @l = split(/\s+/, $_);
		print STDOUT "$1\t$2\t$l[1]\t$l[2]\t" if ($l[0] =~ /([^:]+):(\d+)/);
		for (my $i=3; $i<=$#l; $i += 3) {
			my @pl = @l[$i..$i+2];
			if ($scale) {
				my $sum = 0;
				map {$sum += $_} @pl;
				@pl = map {sprintf "%.9g", $_/$sum} @pl;
			}
			print STDOUT "@pl";
			print STDOUT $i < $#l-2 ? "\t" : "\n";
		}
	}
}

close $stream;

}

sub alleleFreq {

die(qq/
callsetTools.pl alleleFreq [chr number] [position] [sample list] [vcf]

Sample list:
2-column file in sampleID output format where each row lists (1) Sanger ID of sample, (2) group identity.
Assumes a header line is present. Allele frequencies will be calculated for each group.

VCF:
Indexed VCF file.
\n/) unless ((@ARGV && scalar(@ARGV) == 4) || @_);

#my ($chr, $pos, $idlist, $vcf) = @ARGV;
my ($chr, $pos, $idlist, $vcf);
if (@_) {
	($chr, $pos, $idlist, $vcf) = @_;	
} else {
	($chr, $pos, $idlist, $vcf) = @ARGV;
}

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
close $listfh;

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
	my $altfreq = $ndata > 0 ? $altsum / (2*$ndata) : 0;
	print STDOUT "$id\t$altfreq\t$n\t$ndata\n";
}

}

sub multiFreq {
die(qq/
callsetTools.pl multiFreq [-v vcf file name | -V vcf list] [position file] [sample list file] {

Prints ALT allele frequency for each group in the sample list file.

Positon file:
2-column TSV file with fields (1) chromosome and (2) position. Assumes no header.

Sample list file:
2-column file in sampleID output format where each row lists (1) Sanger ID of sample, (2) group identity.
Assumes a header line is present. Allele frequencies will be calculated for each group.

VCF file:
Either the name of a single VCF file using -v or a file listing the VCFs to extract allele frequencies from with -V.
It is assumed that VCFs have names in the format '*_chr1.vcf.gz' if using a list.
\n/) if (!@ARGV || scalar @ARGV < 4);

my $samplefile = pop @ARGV;
die ("Unable to locate sample list file $samplefile\n") if (!-f $samplefile);

my @sites;
my $posfile = pop @ARGV;
open(my $posfh, '<', $posfile) or die("Unable to locate position file $posfile: $!\n");
while (<$posfh>) {
	chomp;
	push @sites, [split(/\s+/, $_)];
}

my %opts=();
getopts('v:V:', \%opts);
my %vcfs;
if ($opts{v}) {
	$vcfs{file} = $opts{v};
	die("Unable to locate VCF file $vcfs{file}\n") if (!-e $vcfs{file});
} elsif ($opts{V}) {
	open(my $vcffh, '<', $opts{V}) or die("Unable to open VCF file list $opts{V}: $!\n");
	while (<$vcffh>) {
		chomp;
		$vcfs{$1} = $_ if ($_ =~ /chr(\d+)\.vcf/);
	}
} else {
	die("VCF file (-v) or VCF file list (-V) required\n");
}
die("No VCFs found\n") if (! %vcfs);

my $header = "CHR\tPOS\tREF\tALT\tAA\tFILTER\tINFO";
my $c = 0;
foreach (@sites) {
	my ($chr, $pos) =@{$_};
	my $chrn = $1 if ($chr =~ /chr(\d+)/);
	my $vcf = exists $vcfs{file} ? $vcfs{file} : $vcfs{$chrn};
	
	my $afout;
	{
		open local(*STDOUT), '>', \$afout;
		alleleFreq($chrn, $pos, $samplefile, $vcf);
	}

	my @tok = split("\n",$afout);
	if ($c == 0) {
		for(my $i=2; $i<=$#tok; $i++) {
			$header .= "\t$1" if ($tok[$i] =~ /^(\S+)/);
		}
		print STDOUT "$header\n";
	}

	my $site = "$chr\t$pos";
	$site .= "\t$1\t$2\t$3\t$4\t$5" if ($tok[0] =~ /ref=(\w)\s+alt=(\w)\s+anc=(\w)\s+(\S+)\s+(\S+)/);
	for (my $i=2; $i <= $#tok; $i++) {
		$site .= "\t$1" if ($tok[$i] =~ /^\S+\s+(\S+)/);
	}
	print STDOUT "$site\n";

	$c++;
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
