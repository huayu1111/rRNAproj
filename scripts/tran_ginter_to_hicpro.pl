#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;

my ($ginterfile,$bedfile,$hicprofile,$help);
GetOptions(
    "ginterfile=s" => \$ginterfile,
    "bedfile=s" => \$bedfile,
    "hicprofile=s" => \$hicprofile,
	"help!" => \$help,
);

my (%giHash,%bedHash);
open(GI,"<$ginterfile") or die "$!\n";
while(<GI>){
    my $line = $_;
    chomp $line;
    my @fieldValues = split /\s+/,$line;
    $giHash{"$fieldValues[0]--$fieldValues[1]--$fieldValues[2]"}{"$fieldValues[3]--$fieldValues[4]--$fieldValues[5]"} = $fieldValues[6];
}
close GI;

open(BED,"<$bedfile") or die "$!\n";
while(<BED>){
    my $line = $_;
    chomp $line;
    my @fieldValues = split /\s+/,$line;
    $bedHash{"$fieldValues[0]--$fieldValues[1]--$fieldValues[2]"} = $fieldValues[3];
}
close BED;

my %outHash;
foreach my $bin1 (keys %giHash){
    foreach my $bin2 (keys %{$giHash{$bin1}}){
       $outHash{$bedHash{$bin1}}{$bedHash{$bin2}} = $giHash{$bin1}{$bin2};
    }
}
close HP;

open(HP,">$hicprofile") or die "$!\n";
foreach my $bin1 (sort {$a <=> $b} keys %outHash){
    foreach my $bin2 (sort {$a <=> $b} keys %{$outHash{$bin1}}){
        print HP "$bin1\t$bin2\t$outHash{$bin1}{$bin2}\n";
    }
}
close HP;