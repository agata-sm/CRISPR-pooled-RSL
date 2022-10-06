#!c:/perl/bin/perl.exe

# Agata Smialowska
# 17 Feb 2022

# script to parse fastq files 
# remove middle part of the fastq read (based on the position)

# usage
#perl fastq_trim_mid_v2.pl --infile in.fastq --pos 24 --len 20
# works with fastq.gz


use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Basename;
	
#commandline parsing for parameters
GetOptions(
	'infile=s'		=>	\(my $path2infile),		#path/to/input/fastq_file
	'pos=s'		=>	\(my $pos),		#position of the 1st base to be removed
	'len=s'		=>	\(my $seqlen),		#length of the sequence to be removed
) or die "Error in command line arguments";



my($infile_name, $dirs, $suffix) = fileparse($path2infile);
print "processing file $infile_name\n";

my $gz;
if ($path2infile=~m/\S+\.fastq$/){
	$gz="no";
}
elsif ($path2infile=~m/\S+\.fastq\.gz$/){
	$gz="yes";
}


$infile_name=~m/(^\S+)\.fastq/;
my $library_name=$1;
print "$library_name\n";


my $out_file_trimmed="$library_name\.trimmed.fastq";

my $counter=0;
my $counter_trimmed=0;


open (INFILE, "<", $path2infile) or die "Cannot open file $path2infile: $!"; 
open (OUTFILE_T, ">", $out_file_trimmed) or die "Cannot open logfile file $out_file_trimmed: $!"; 

if($gz eq qw /no/){
	open (INFILE, "<", $path2infile) or die "Cannot open file $path2infile: $!"; 
}elsif($gz eq qw /yes/){
	open (INFILE, "gzip -dc < $path2infile |") or die "Cannot open file $path2infile: $!"; 
}




#read in fastq file in chunks 4 lines each (1 record)
while ((my @lines = map $_ = <INFILE>, 1..4) [0]) {
	chomp($lines[0],$lines[1],$lines[2],$lines[3]);

	my @seq_line=split //, $lines[1];
	splice @seq_line, $pos-1, $seqlen; # positions $pos to $pos+$seqlen (0-based indices)
	my $seqline_trimmed=join '',@seq_line;

	my @qvals=split //, $lines[3];
	splice @qvals, @-, $seqlen;
	my $qvals_trimmed=join '',@qvals;


		print OUTFILE_T "$lines[0]\n$seqline_trimmed\n$lines[2]\n$qvals_trimmed\n";

		$counter_trimmed++;
		$counter++;
	}
	

print "reads processed: $counter\ntrimmed: $counter_trimmed\n\n";

