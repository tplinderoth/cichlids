#/usr/bin/perl

while ($arg = shift) {
    if ($arg =~ s/^-//) { $opt{$arg} = 1 ; }
    else { last ; }
}

$chr = $arg ;
($pos = shift) || die "usage: perl looksite.pl [-ref] <chr> <pos> [<field>] where field is simple_id, genus, species, clade\n    need to set HEADER, METADATA, VCFDIR\n" ;
$feat = shift ;
($chr >= 1 && $chr <= 23) || die "error: chr $chr must be between 1 and 23\n" ;

open (H, "<HEADER") || die "error: can't open HEADER\n" ;
$_ = <H> ; chomp ; ($c, $p, $id, $ref, $alt, $qual, $filt, $info, $format, @pid) = split /\t/ ;
close H ;

open (M, "<METADATA") || die "error: can't open METADATA\n" ;
$_ = <M> ; chomp ; ($p, @cols) = split /\t/ ; 
pop @cols ; 			# hack - for some reason there is a final empty column
for ($i = 0 ; $i < @cols ; ++$i) { $col{$cols[$i]} = $i ; }
($ic = $col{"simple_id"}) ;
($gc = $col{"genus"})      || die "error: can't find genus column\n" ;
($sc = $col{"species"})    || die "error: can't find species column\n" ;
while (<M>) {
    chomp ; ($p, @fields) = split /\t/ ;
    $fields[$sc] = $fields[$gc] . " " . $fields[$sc] ;
    $meta{$p} = [@fields] ;
}
close M ;

 # defined $meta{$pid} || die "can't find $pid in meta\n" ;
 # defined $col{$feat} || die "can't find column for $feat\n" ; 
 # print "$feat for $pid is $meta{$pid}[$col{$feat}]\n" ;

for $p (@pid) { defined $meta{$p} || die "error: can't find $p in metadata\n" ; }

if ($feat && !defined $col{$feat}) { die "error: feature $feat is not defined\n" ; }
$fc = $col{$feat} ;

open (V, "bcftools view -H VCFDIR/malawi_cichlid_v1_chr$chr.vcf.gz chr$chr:$pos |") || die "error: can't access VCF\n" ;
while (<V>) {
    chomp ; ($c, $p, $id, $ref, $alt, $qual, $filt, $info, $format, @g) = split /\t/ ;
    print "$c $p $ref $alt $info\n" ;
    for ($i = 0 ; $i < @pid ; ++$i) {
	if ($feat) { $k = $meta{$pid[$i]}[$fc] ; } else { $k = "all" ; }
	$g[$i] =~ /(.)\/(.)/ ; $g1 = $1 ; $g2 = $2 ;
	++$t{$k} if ($g1 ne ".") ; ++$t{$k} if ($g2 ne ".") ;
	if ($opt{"ref"}) {
	    ++$c{$k} if ($g1 == 0) ; ++$c{$k} if ($g2 == 0) ;
	} else {
	    ++$c{$k} if ($g1 == 1) ; ++$c{$k} if ($g2 == 1) ;
	}
    }
    foreach $k (sort keys %c) {
	print ("$k\t$c{$k}\t$t{$k}\n") ;
    }
}
