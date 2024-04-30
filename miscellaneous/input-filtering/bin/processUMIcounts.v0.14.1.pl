#!c:/perl/bin/perl.exe

# author: Agata Smialowska

# NBIS proj 5351

# TESTED for perl (v5.18.4)


# script to filter and normalise UMI counts
# consolidate barcodes per guide (barcoded screens) i.e. count number of unique barcodes per guide


#############################
# input files

# count table - from CrisprCounter.jar
# RSL.guide,guide.set,smpl1,smpl2,smpl3,smpl4
# GM2A_4_AGATTTTTCA,GM2A_4,0,1,0,0
# GM2A_4_TATGGATTAA,GM2A_4,0,0,0,1
# GM2A_4_TGTAAATGAC,GM2A_4,1,0,0,0


# library confing file
# Guide,Sequence,Gene
# A1BG_1,CATCTTCTTTCACCTGAACG,A1BG
# A1BG_2,CTCCGGGGAGAACTCCGGCG,A1BG
# A1BG_3,TCTCCATGGTGCATCAGCAC,A1BG
# A1BG_4,TGGAAGTCCACTCCACTCAG,A1BG


# prefiltered input (output of filter_RSL_input.v0.10.pl)
# RSL_guide,guide_set,InputFilt
# GM2A_4_CCAATGTATA,GM2A_4,15
# GM2A_4_GGACAAGAAA,GM2A_4,10
# GM2A_4_AACCCACCGT,GM2A_4,9



#############################
# filtering

# filtering out artefacts ("ghosts") - RSL-sgRNA combinations with very few reads
# based on the presence in the prefiltered master input file (output of filter_RSL_input.v0.10.pl)
# additional filter for RowSums, set manually (to remove observations with very low count, e.g. 1 across all replicates)



#############################
# normalisation

# norm factor
# average of all the column sums (AvColSum); calculate the normalization factor AvColSum/ColSum for each column;
# multiply each column with its normalization factor 


#############################
# OUTPUT

# diagnostic
# output of frequencies of UMI read counts (not normalised, after all filtering steps) in each sample: for histogram plotting

# data
# count tables 
# unique RSL / UMI counts per guide (counting after filtering and normalisation)
# read counts (unfiltered, not normalised) per guide

# the output is saved in separate directories for read counts and RSL/UMI counts


#############################
# new in version 0.12

# fix the read counting - is now done before the input based filter
# input based filter is optional 

# new in version 0.13

# output is simplified, all files saved in one directory given as the command line argument
# first column in the output should be guide not RSL.guide!

# new in v 0.14
# use csv library definition

# v0.14.1 adds argument with outfiles prefix > necessary for input filtering pipeline

use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use List::Util qw( sum );
use File::Basename;
use File::Path qw(make_path);
#use Data::Dumper qw(Dumper);

use experimental 'smartmatch';
no warnings 'experimental::smartmatch'; ### Careful with possible changes in implementation in perl > 5.18

my $script_name="processUMIcounts.v0.14.1.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/SAMPLE.UMIcounts.csv\n";
	print "--outdir: /path/to/outdir\n";
	print "--prefix: prefix for outfiles\n";
	print "--input_lib: /path/to/processed_input_library NAME.filtered.csv (optional; if not given, no input filtering is performed)\n";	#Input_Brunello_mod.filtered.csv
	print "--input_lib_design: /path/to/file with library design Library.csv; comma delimited: Guide - Sequence - Gene\n";	#Brunello_Library_USE_THIS_ONLY.csv
	print "--filter: select filter for removing the 0-inflated part of data; format: ´CO=n´ where n is the cut-off value for rowSums;\nsgRNA-RSL combinations with (library size normalised) rowSum counts larger or equal to cut-off are included;\nif no rowSum filter is desired please use CO=0;\n";
	print "-h prints this message\n";
}


else{
	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outdir=s'		=>	\(my $outdir),
		'pref=s'		=>	\(my $pref),
		'input_lib_design=s'		=>	\(my $input_lib_design),
		'input_lib:s'		=>	\(my $input_lib),			
		'filter=s'		=>	\(my $filter)
	) or die "Error in command line arguments";


	unless (-d $outdir) {
    	make_path($outdir);
	}


	my @suffixlist=(".tsv",".csv",".txt");
	#my $basename = basename($infile,@suffixlist);

	my $basename=$pref;
	
	print "$basename\n\n";

	my $outfile_UMI_perguide_fname="$basename\.RSL.perguide.tsv";
	my $outfile_reads_all_perguide_fname="$basename\.nonfilt_reads.perguide.tsv";
	my $outfile_freqs_raw_fname="$basename\.frequencies\.filt_reads.tsv";


	my $log_fname="$basename\.readme.log"; #readme with cmd and stdout echo

	my $outfile_UMI_perguide="$outdir\/$outfile_UMI_perguide_fname";
	my $outfile_reads_all_perguide="$outdir\/$outfile_reads_all_perguide_fname";
	my $outfile_freqs_raw="$outdir\/$outfile_freqs_raw_fname";


	my $log="$outdir\/$log_fname";

	my $timenow=localtime();
	open (LOG, ">", $log) or die "Cannot open logfile file $log: $!"; 
	print LOG "This is log file for analysis using $script_name on $timenow\n";

	open (STDOUT, "| tee $log") or die "Teeing off: $!\n";

	print "This is log file for analysis using $script_name on $timenow\n";


	print "Project name is $basename\n\n";
	print "parameters: $parameters\n";
	print "\n\nAll output files $basename `.*` in directory RSL are of reads filtered based on RowSum filter $filter \n";

	if(length $input_lib){
		print "and presence in the prefiltered Input library $input_lib";
	}

	

	print "\n\noutput files are\n";

	print "\nDiagnostic:\n";

	print "$outfile_freqs_raw_fname - frequencies of non-normalised read counts per RSL-guide after filtering (for presence in input and rowSums, if applicable)\n";

	print "\nOutput:\n";

	print "Output directory is $outdir\n";

	print "$outfile_reads_all_perguide_fname - read counts summarised per guide, no input filtering applied, no post-filtering applied, no RSL evaluation applied\n";
	print "$outfile_UMI_perguide_fname - RSL counts per guide after filtering\n";



	print "\nLog:\n";
	print "$log_fname - log file which includes this text\n\n";



	#check filter argument and set the cutoff value
	my $filt_cutoff;
	if ($filter =~ m/CO=(\d+\.?\d*)$/ ){
		$filt_cutoff=$1;
	}else{
		print "Selected cutoff CO is not a number. The correct format is CO=n where n is a number (floating point accepted).\n";
		print "No filter used\n";
		$filt_cutoff=0;
	}

	print "Selected filtering strategy is ";

	if (length $input_lib){
		print "prefiltering based on prefiltered Input library $input_lib and ";
		}
	print "RowSum filter $filter and the cutoff value is $filt_cutoff\n\n";



	#process the reads from processed filtered input library
	my %input_RSL_guide;
	if (length $input_lib){

		open (INFILE_INPUT, "<", $input_lib) or die "Cannot open file with processed input library $input_lib: $!";
		while (<INFILE_INPUT>){
			chomp $_;
			my @line=split /,/;	

			unless ($_ =~m/RSL\.guide/){

				my $RSL_guide=shift @line;
				my $guide_set=shift @line;

				$input_RSL_guide{$RSL_guide}=$line[0];

			}

		}	

		close(INFILE_INPUT);
	}


	# gene names from library design file
	my %gene_guide;
	open (INFILE_INPUT_LIBDES, "<", $input_lib_design) or die "Cannot open file with processed input library $input_lib_design: $!";
	while (<INFILE_INPUT_LIBDES>){
	    # change line endings to \n from \r\n #### OBS! added this here, absent in the main pipeline
 		$_=~s/[\r\n]+//;

		chomp $_;

		unless ($_ =~m/Guide/){ #this file is header-less but just in case
			my @line=split /,/;	

			my $guide=$line[0];
			my $gene=$line[2];

			$gene_guide{$guide}=$gene;

		}

	}	
	close(INFILE_INPUT_LIBDES);


	########################################
	### first pass over the count table - to calculate normalisation factors and get reads summarised per guide

	my $timenow2=localtime();
	print "\n\n$timenow2: Calculate normalisation factors\n";


	#calc norm factors
	my @col_sums;
	my $number_of_samples;

	#raw read counts for outputing per guide in $outfile_reads_all_perguide
	my %guide_raw_read_counts; #HoA; the keys are $RSL_guide combination
	
	#header
	my $header_per_guide_raw;
	my $header_smpls;



	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 

	while(<INFILE>){

		chomp $_;
		my @line=split /,/;


		unless ($_ =~m/RSL\.guide/){

			my $RSL_guide=shift @line;
			my $guide_set=shift @line;

			$number_of_samples=scalar(@line);

			# col sums calculation AFTER filtering for presence in input, if defined
			if (length $input_lib){
				if ( exists($input_RSL_guide{$RSL_guide}) ){
					foreach my $i (0.. $#line) {
						$col_sums[$i] += $line[$i];	
					}
				}

			}else{ #input library not defined
				foreach my $i (0.. $#line) {
					$col_sums[$i] += $line[$i];	
				}
			}

			###############
			# get reads summarised per guide

			if(exists($guide_raw_read_counts{$guide_set})){
				my @counts_guideset=@{$guide_raw_read_counts{$guide_set}};
					
				my @updated_counts_guideset;
				foreach my $i (0 .. $#counts_guideset){

		 	 		$updated_counts_guideset[$i]=$counts_guideset[$i]+$line[$i];
		 	 	}
		 	 	$guide_raw_read_counts{$guide_set}=[@updated_counts_guideset];

			}else{
				$guide_raw_read_counts{$guide_set}=[@line];
			}


		}else{ #get header
			my $header1=shift @line;
			my $header2=shift @line;
			my @UMI_counts=@line;

			$header_smpls=join ",", @UMI_counts;

			my $line_out="$header1,$header2,$header_smpls";

			$header_smpls=~s/,/\t/g;
			$header_per_guide_raw="sgRNA\tGene\t$header_smpls";


		}
	}


	close(INFILE);

	if (length $input_lib){
		print "\nNormalisation factors are calculated AFTER input based filtering\n";
	}

	print "Col sums are @col_sums\n";

	my $tot_sum=sum(@col_sums);


	my $AvgColSum=$tot_sum/$number_of_samples;

	print "Total sum is $tot_sum\nAverage sum is $AvgColSum\tNumber of samples is $number_of_samples\n";
	my @norm_factors;
	
	foreach my $i (0 .. $#col_sums){

		$norm_factors[$i]=$AvgColSum/$col_sums[$i];

	}
	print "Normalisation factors are\n@norm_factors\n";



	################################
	#### output non filtered reads summarised per guide

	#output the file with raw read counts
	open (OUTFILE_RAW, ">", "$outfile_reads_all_perguide") or die "Cannot open output file $outfile_reads_all_perguide: $!";
	print OUTFILE_RAW "$header_per_guide_raw\n";

	for my $guide_set ( sort keys %guide_raw_read_counts ){
		my $gene=$gene_guide{$guide_set};
		
		my @total_counts_per_guide=@{$guide_raw_read_counts{$guide_set}};
		
		my $tot_counts_for_file=join "\t",@total_counts_per_guide;
		print OUTFILE_RAW "$guide_set\t$gene\t$tot_counts_for_file\n";

	}


	close(OUTFILE_RAW);

	undef %guide_raw_read_counts;





	################################
	### normalise, filter and summarise by guide


sub normalise_counts{

	my ($line_ref, $norm_factors_ref) = @_;

	my @line_sub=@$line_ref;
	my @norm_factors_sub=@$norm_factors_ref;

	my @norm_counts;

	foreach my $i (0 .. $#line_sub){

		#$norm_counts[$i]=$line[$i]*$norm_factors[$i];
		$norm_counts[$i]=sprintf ("%.3f",$line_sub[$i]*$norm_factors_sub[$i]);

	}

	return @norm_counts;
}


# sub count_UMIs{

# 	my ($guide_UMI_counts_ref, $norm_counts_ref, $guide_set)=@_;

# 	my @norm_counts_sub=@$norm_counts_ref;
# 	my %guide_UMI_counts_sub=%$guide_UMI_counts_ref;

# 	my @UMI_counts= map {if ($_ > 0) { 1 } else { 0 } } @norm_counts_sub;

# 	if(exists( $guide_UMI_counts_sub{$guide_set})){
# 		my @counts_guideset=@{$guide_UMI_counts_sub{$guide_set}};
# 		my @updated_counts_guideset;
		
# 		foreach my $i (0 .. $#counts_guideset){
# 			$updated_counts_guideset[$i]=$counts_guideset[$i]+$UMI_counts[$i];
# 		}
# 		#update the original hash not the local subroutine copy
# 		$ {$guide_UMI_counts_ref}{$guide_set}=[@updated_counts_guideset];

# 	}else{
# 		$ {$guide_UMI_counts_ref}{$guide_set}=[@norm_counts_sub];
# 	}

# 	%{$guide_UMI_counts_sub} = ();

# }					


sub count_UMIs{

	my ($guide_UMI_counts_ref, $norm_counts_ref, $guide_set)=@_;

	my @norm_counts_sub=@$norm_counts_ref;

	my @UMI_counts= map {if ($_ > 0) { 1 } else { 0 } } @norm_counts_sub;

	if(exists( ${$guide_UMI_counts_ref}{$guide_set})){
		my @counts_guideset=@{${$guide_UMI_counts_ref}{$guide_set}};
		my @updated_counts_guideset;
		
		foreach my $i (0 .. $#counts_guideset){
			$updated_counts_guideset[$i]=$counts_guideset[$i]+$UMI_counts[$i];
		}
		$ {$guide_UMI_counts_ref}{$guide_set}=[@updated_counts_guideset];

	}else{
		$ {$guide_UMI_counts_ref}{$guide_set}=[@UMI_counts];
	}


}	

	my $timenow3=localtime();
	print "\n\n$timenow3: Filter and summarise by sgRNA\n";

	# for summarising UMI / RSL counts per sgRNA
	my %guide_UMI_counts; #HoA
	my $header_per_guide;


	# for read count frequencies

	# get the values for read counts to be able to initialise the hash with arrays containing elements "0" of length equal to number of samples
	# this is array with the values of counts, which are not incremented by 1 from [0, max] especially for higher counts; holds the "seen" values
	my @values_raw_readcounts;

	#frequencies
	my %freq_raw_reads; #HoH where the keys of the outer hash are the samples and the inner hash is: RAW read count - frequency


	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 

	while(<INFILE>){

		chomp $_;
		my @line=split /,/;

		my $header;

		if ($_ =~m/RSL\.guide/){ #header line

			#samples
			my $header1=shift @line;
			my $header2=shift @line;
			my @UMI_counts=@line;

			my $header_smpls=join ",", @UMI_counts;

			my $line_out="$header1,$header2,$header_smpls";

			$header_smpls=~s/,/\t/g;
			$header_per_guide="sgRNA\tGene\t$header_smpls";

			print "Sample names are\n$header_smpls\n";


		}else{ #line with data

			my $RSL_guide=shift @line;
			my $guide_set=shift @line;

			#filtering for presence in input library, if defined
			if (length $input_lib){
				if ( exists($input_RSL_guide{$RSL_guide}) ){

					my @norm_counts=normalise_counts(\@line,\@norm_factors);
					my $norm_sum=sum(@norm_counts);
				
					if($norm_sum >= $filt_cutoff){

						count_UMIs(\%guide_UMI_counts,\@norm_counts,$guide_set);

						foreach my $i (0.. $#line) {

							unless ($line[$i] ~~ @values_raw_readcounts){
								push @values_raw_readcounts,$line[$i];
							}
					
							if(exists($freq_raw_reads{$i}{$line[$i]} )) {
								$freq_raw_reads{$i}{$line[$i]}++;
							}else{
								$freq_raw_reads{$i}{$line[$i]}=1;
							}
						}					
					}
				}
			}else{ #input library not defined
				
				my @norm_counts=normalise_counts(\@line,\@norm_factors);
				my $norm_sum=sum(@norm_counts);
				
				if($norm_sum >= $filt_cutoff){

					count_UMIs(\%guide_UMI_counts,\@norm_counts,$guide_set);
		
					foreach my $i (0.. $#line) {

						unless ($line[$i] ~~ @values_raw_readcounts){
							push @values_raw_readcounts,$line[$i];
						}
					
						if(exists($freq_raw_reads{$i}{$line[$i]} )) {
							$freq_raw_reads{$i}{$line[$i]}++;
						}else{
							$freq_raw_reads{$i}{$line[$i]}=1;
						}
					}
				}			
			}
		}
	}#proc infile
	close(INFILE);

	open (OUTFILE_GUIDE, ">", "$outfile_UMI_perguide") or die "Cannot open output file $outfile_UMI_perguide: $!";
	print OUTFILE_GUIDE "$header_per_guide\n";

	#consolidation of counts per guide set
	for my $guide_set ( sort keys %guide_UMI_counts ){

		my $gene= $gene_guide{$guide_set};

		my @total_counts_per_guide=@{$guide_UMI_counts{$guide_set}};
		my $tot_counts_for_file=join "\t",@total_counts_per_guide;
		print OUTFILE_GUIDE "$guide_set\t$gene\t$tot_counts_for_file\n";


	}

	close(OUTFILE_GUIDE);


	################################
	### output frequencies of reads per sgRNA-RSL combination after filtering


	# build a printable structure for the frequencies data: HoA where keys are numbers of UMI-sgRNA reads from @values_raw_readcounts 
	# and values are frequencies in each sample
	my %proc_freq_raw_reads;
	#init the arrays for each seen readcount value
	my @mock_counts=("0") x $number_of_samples;

	for my $readcount_seen (@values_raw_readcounts){
		$proc_freq_raw_reads{$readcount_seen}=[@mock_counts];
	}

	for my $sample_id ( sort {$a <=> $b} keys %freq_raw_reads ){
		foreach my $reads_number (sort {$a <=> $b} keys %{$freq_raw_reads{$sample_id}}) {

    		${$proc_freq_raw_reads{$reads_number}}[$sample_id]=$freq_raw_reads{$sample_id}{$reads_number};

    	}

	}

	#frequencies file
	open (OUTFILE_FREQ_RAW, ">", "$outfile_freqs_raw") or die "Cannot open output file $outfile_freqs_raw: $!";
	my $header_freqs="readcount\t$header_smpls";
	print OUTFILE_FREQ_RAW "$header_freqs\n";

	for my $readcount (sort {$a <=> $b} keys %proc_freq_raw_reads){
		
		my @frequencies_in_libraries=@{ $proc_freq_raw_reads{$readcount} };

		my $frequencies_to_print=join "\t",@frequencies_in_libraries;
		print OUTFILE_FREQ_RAW "$readcount\t$frequencies_to_print\n";

	}


	close (OUTFILE_FREQ_RAW);

	undef %proc_freq_raw_reads;
	undef %freq_raw_reads;



}

my $timenow4=localtime();
print "*** $timenow4: Finished! ***\n";
exit;


