#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;

my ($eigenfile,$outfile,$help);
GetOptions(
    "eigenfile=s" => \$eigenfile,
    "outfile=s" => \$outfile,
	"help!" => \$help,
);

open(EI,"<$eigenfile") or die "$!\n";
open(OUT,">$outfile") or die "$!\n";
print OUT "track name=\"CH12LX PC1\"  yLineMark=\"0.0\" alwaysZero=on maxHeightPixels=100:75:11 visibility=full color=255,128,0 altColor=51,51,255 viewLimits=-1:1 autoScale=on type=bedGraph";

my $len = 98319150;
my $itererNum = 1;
while(<EI>){
    my $line = $_;
    chomp $line;
    if($itererNum < 197){
        my $start = ($itererNum-1)*500000+1;
        my $end = $itererNum*500000;
        print OUT "chr16\t$start\t$end\t$line\n";
    }else{
        my $start = ($itererNum-1)*500000+1;
        my $end = $len;
        print OUT "chr16\t$start\t$end\t$line\n";
    }
    $itererNum = $itererNum + 1;
}
close EI;
close OUT;
