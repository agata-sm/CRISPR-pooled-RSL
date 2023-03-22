#!c:/perl/bin/perl.exe

#script to parse stdout of UMI summarising script CrisprCounter.jar


use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use List::Util 'first';  

my $script_name="parse_CrisprCounter.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/counter.stdout.txt\n";
	print "--outfile: /path/to/outfile counter.stdout.summary.txt\n";
}else{

	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outfile=s'		=>	\(my $outfile)
	) or die "Error in command line arguments";



	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 

	undef $/;

	my @chunks = split("Reading guide fastq file ", <INFILE>);

	close INFILE;

	my $config=shift @chunks;

	open (OUTFILE, ">","$outfile") or die "Cannot open output file $outfile: $!"; 

	my $header="file_fastq\tsgRNA_file\tsgRNA_total\treads_assigned_file\treads_total\tfraction_reads_assigned";

	print OUTFILE "$header\n";

	foreach my $record (@chunks){

		my @lines=split("\n",$record);

		my $file_fastq = first { /.fastq/ } @lines;
		my $gRNA_line= first{ /guide RNAs observed/} @lines;
		my $mapped_line= first {/mapped to guideRNA/} @lines;
		my $reads_line= first {/reads read in total/} @lines;

		#remove leading path to fastq file
		$file_fastq =~ s{^.*/}{}; 
		
		my @file_stats;

		if($gRNA_line =~m/^(\d+) out of (\d+) guide RNAs observed/){
			
			push @file_stats, $1;
			push @file_stats, $2;

		}

		if($mapped_line =~m/^(\d+) mapped to guideRNA/){
			push @file_stats, $1;
		}
		if($reads_line =~m/^(\d+) reads read in total/){
			push @file_stats, $1;
		}

		my $perc_reads_mapped=sprintf ("%.3f",$file_stats[2]/$file_stats[3]);
		push @file_stats, $perc_reads_mapped;
		
		my $file_stats_line=join "\t", @file_stats;

		print OUTFILE "$file_fastq\t$file_stats_line\n";

	}
close OUTFILE;

}

exit;

