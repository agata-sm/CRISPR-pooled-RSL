#!c:/perl/bin/perl.exe

# author: Agata Smialowska

# NBIS proj 5351

# TESTED for perl (v5.18.4)


# script to calculate log2FC of UMI counts and rank them by the log2FC; ties are broken by average UMI count
# used for barcoded CRISPR screens

# input file:
# NAME.RSL.perguide.tsv
# produced by processUMIcounts.pl

# output file
# NAME.RSL.log2FCrank.tsv
# gene - log2FC.condition - pos.Rank - neg.Rank



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use List::Util qw( sum max );
use File::Basename;
use File::Path qw(make_path);
#use Data::Dumper qw(Dumper);
use Sort::Rank qw(rank_sort rank_group);


my $script_name="rank_log2FC.v0.1.1.pl";



if ($ARGV[0]eq qw '-h'){
	print "this script ranks the guides based on log2FC in treatment vs reference\n";
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/RSL.perguide.tsv\n";
	print "--outfile: /path/to/outfile\n";
	print "--ref: name of the reference sample (colname in the input file)\n";
	print "--treat: name of the treatment sample (colname in the input file)\n";
	print "-h prints this message\n";
}


else{
	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outfile=s'		=>	\(my $outfile),
		'ref=s'		=>	\(my $ref_smpl),
		'treat=s'		=>	\(my $treat_smpl)			
	) or die "Error in command line arguments";



	my @suffixlist=(".tsv",".csv",".txt");
	my $basename = basename($infile,@suffixlist);

	print "$basename\n\n";

	# pseudocount to avoid division by 0
	my $pseudocnt=0.01;

	# Function for log2 calculator
	sub log2 
	{
   	 	my $n = shift;
      
   		# using pre-defined log function
    	return log($n) / log(2);
	}

	#logfile to record the contrast for log2FC
	my $outdirname  = dirname($outfile);

	print "Outdir is $outdirname\n";

	my $logfile="$outdirname\/$basename\.log";

	print "Log file is $logfile\n";

	open (LOG, ">","$logfile") or die "Cannot open log file $logfile: $!"; 
	print LOG "Ranks saved in file\n$outfile\n\nare based on contrast\n\n$treat_smpl vs. $ref_smpl\n\ndata count table\n$infile\n\n";
	print LOG "Pseudocount $pseudocnt is added to the counts to avoin division by 0\n";
	close(LOG);


	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 


	#get the indices of the columns with samples of interest
	my @treat_idx;
	my @ref_idx;
	my $treat_idx;
	my $ref_idx;

	my %guide_logFC;
	my %guide_aveCount;

	my %hash_foo;
	my @scores_foo;

	while(<INFILE>){

		chomp $_;
		my @line=split /\t/;

		
		if ($_ =~m/sgRNA/){
			@treat_idx=grep { $line[$_] =~ /$treat_smpl/ } 0..$#line; 
			@ref_idx=grep { $line[$_] =~ /$ref_smpl/ } 0..$#line; 
		
			$treat_idx=$treat_idx[0];
			$ref_idx=$ref_idx[0];


		}
		else{ #data

			my $cnt_treat=$line[$treat_idx[0]] + $pseudocnt;
			my $cnt_ctrl=$line[$ref_idx[0]] + $pseudocnt;

			my $log2FC=log2( $cnt_treat/$cnt_ctrl );
			my $averageCount=($line[$treat_idx[0]] + $line[$ref_idx[0]]) / 2;

			$guide_logFC{$line[0]}=$log2FC;

			$guide_aveCount{$line[0]}=$averageCount;

			my %hash_foo;
			$hash_foo{'score'}=$log2FC;

			my $rec={};
			$rec->{'score'} = $log2FC;
			$rec->{'sgRNA'} = $line[0];

			push @scores_foo, $rec;

		}
	}

	close(INFILE);


	my @sorted = rank_group(\@scores_foo);


#resolve this structure

# Given an array reference and optional score extraction subroutine return 
# an array containing the elements of the input array arranged in rank order. 
# Each element of the returned array is a reference to a three element array 
# containing the rank of this element, a flag that indicates whether this rank is shared with other elements and the corresponding value from the input array.

# save the ranked results to a hash where keys are sgRNAs >> this is because the reverse rank is also output in the same file
# calculate the rev rank when printing the output: calculate reverse of the ranks i.e. list where lowest -logFC get highest ranks
# rev rank = max rank - rank +1



my @ranks;

my %ranks_per_sgRNA;

for (@sorted){

	my $sorted_lvl1_ref=$_;

	my @sorted_lvl1=@$sorted_lvl1_ref;

	my $number_of_hashes=scalar(@sorted_lvl1)-1;

	my $rank=$sorted_lvl1[0];

	push @ranks, $rank; 

	#print the contents of each hash to the file
	foreach my $i (1..$number_of_hashes){
			
		my $hash_lvl2_ref=$sorted_lvl1[$i];
		my %hash_lvl2=%$hash_lvl2_ref;

		my $score_2prnt=$hash_lvl2{'score'};;
		my $sgRNA_2prnt=$hash_lvl2{'sgRNA'};

		#create a new hash with ranks to later calculate reverse of the ranks i.e. list where lowest -logFC get highest ranks
		my $sgRNA_score="$sgRNA_2prnt::$score_2prnt";
		$ranks_per_sgRNA{$sgRNA_score}=$rank;

	}

}

my $max_rank=max(@ranks);

print "Max rank is $max_rank\n";

open (OUTFILE, ">", $outfile) or die "Cannot open outfile $outfile: $!\n";
my $header="sgRNA\tlogFC\trank_poslFC\trank_neglFC";
print OUTFILE "$header\n";

foreach my $sgRNA_score (sort keys %ranks_per_sgRNA) {

	my ($sgRNA,$score)=split/::/,$sgRNA_score;
	my $rank_poslFC=$ranks_per_sgRNA{$sgRNA_score};
	my $rank_neglFC=$max_rank-$rank_poslFC+1;

	print OUTFILE "$sgRNA\t$score\t$rank_poslFC\t$rank_neglFC\n";

}

close OUTFILE;


}
exit;






