#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$help);

GetOptions(
    "inputdir=s" => \$inputdir,
	"help!" => \$help,
);

my @bedgraphfiles = `find $inputdir -name "*bw"`;

foreach my $bedgraphfile (@bedgraphfiles){
    chomp $bedgraphfile;
    my $bwfile = $bedgraphfile;
    $bwfile =~ s/bw$/bedgraph/;
    system("/public/ZhangJinLab/bigWigToBedGraph $bedgraphfile $bwfile")==0 or die "$!\n";
}

