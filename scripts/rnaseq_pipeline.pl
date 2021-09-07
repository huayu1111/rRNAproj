#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$indexdir,$hisat2dir,$picarddir,$stringtiedir,$gfffile,$threads,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
	"indexdir=s" => \$indexdir,
	"hisat2dir=s" => \$hisat2dir,
	"picarddir=s" => \$picarddir,
	"stringtiedir=s" => \$stringtiedir,
	"gfffile=s" => \$gfffile,
	"threads=s" => \$threads,
	"help!" => \$help,
);

my @sample_p1s = `find $inputdir -name "*_R1.fq.gz"`;
print join("\n",@sample_p1s)."\n";
foreach my $sample_p1 (@sample_p1s){
	chomp $sample_p1;
	$sample_p1 =~ /.*\/(.*)\_R1.fq.gz/;
	my $sample_id = $1;
	my $sample_p2 = $sample_p1;
	$sample_p2 =~ s/R1/R2/;
	if(!-e "$outputdir/$stringtiedir/$sample_id"){
		mkpath("$outputdir/$stringtiedir/$sample_id",0644);
		if($@){
			print "Make path $outputdir/$stringtiedir/$sample_id failed:\n";
			exit(1);
		}
	}
	if(!-e "$outputdir/$hisat2dir/$sample_id"){
		mkpath("$outputdir/$hisat2dir/$sample_id",0644);
		if($@){
			print "Make path $outputdir/$hisat2dir/$sample_id failed:\n$@";
			exit(1);
		}
	}
	open(SH,">$outputdir/$hisat2dir/$sample_id/${sample_id}_expcal.sh") or die "$!\n";
	if(!-e "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
		print SH "hisat2 -x $indexdir -p $threads --dta --rg-id $sample_id --rg SM:$sample_id -1 $sample_p1 -2 $sample_p2  --summary-file $outputdir/$hisat2dir/$sample_id/mapping_summary.txt -S $outputdir/$hisat2dir/$sample_id/accepted_hits.sam\n";
	}
	if(!-e "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
		print SH "grep -v -E -w 'NH:i:2|NH:i:3|NH:i:4|NH:i:5|NH:i:6|NH:i:7|NH:i:8|NH:i:9|NH:i:10|NH:i:11|NH:i:12|NH:i:13|NH:i:14|NH:i:15|NH:i:16|NH:i:17|NH:i:18|NH:i:19|NH:i:20' $outputdir/$hisat2dir/$sample_id/accepted_hits.sam > $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sam\n";
        print SH "samtools view -bS $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sam -o $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.bam\n";
	}
	if(!-e "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
		print SH "samtools sort -@ $threads -o $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.bam $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.bam\n";
	}
	if(!-e "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
		print SH "java -Xmx15g -jar $picarddir/picard.jar MarkDuplicates I=$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.bam O=$outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam METRICS_FILE=$outputdir/$hisat2dir/$sample_id/${sample_id}.metricsFile VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate\n";
        print SH "samtools index $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam\n";
        print SH "rm $outputdir/$hisat2dir/$sample_id/accepted_hits_NHi1.sam\n";
        print SH "rm $outputdir/$hisat2dir/$sample_id/accepted_hits.sam\n";
	}

	if(!-e "$outputdir/$stringtiedir/$sample_id/transcripts.gtf" || -z "$outputdir/$stringtiedir/$sample_id/transcripts.gtf"){
		print SH "stringtie -p $threads -e -B -G $gfffile -A $outputdir/$stringtiedir/$sample_id/gene_abund.tab -o $outputdir/$stringtiedir/$sample_id/transcripts.gtf $outputdir/$hisat2dir/$sample_id/accepted_hits.sorted.bam\n";
	}
	close SH;
    
	close SH;
    my $taskNum =`ps -aux | grep stringtie | wc -l`; 
    while($taskNum > 10){
        print "The num of task remaining $taskNum\n";
        sleep 30;
        print `date`;
        $taskNum = `ps -aux | grep stringtie | wc -l`;
    }
	my $out = system("sh $outputdir/$hisat2dir/$sample_id/${sample_id}_expcal.sh 1>>$outputdir/$hisat2dir/$sample_id/std.log 2>>$outputdir/$hisat2dir/$sample_id/error.log &");
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}							