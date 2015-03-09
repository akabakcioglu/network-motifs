#/usr/bin/perl -w
use strict;

my $f = $ARGV[0];
open(FP,"$f.nodeids") or die "Can't open $f.nodeids\n";
my %id=();
while (<FP>){
    chomp;
    my @a = split/\t/;
    $id{$a[0]}=$a[1];
}
close(FP);

open(FP,"$f") or die "Can't open $f\n";
my $i=0;
while (<FP>){
    $i += 1;
    chomp;
    my @a = split/\t/;
    for my $j (1..@a) {
	print join("\t",$id{$j},"+",$id{$i})."\n" if ($a[$j-1]==1);
	print join("\t",$id{$j},"-",$id{$i})."\n" if ($a[$j-1]==-1);
    }
}
close(FP)
