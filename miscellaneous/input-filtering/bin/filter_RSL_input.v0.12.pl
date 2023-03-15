#!c:/perl/bin/perl.exe

# author: Agata Smialowska

# NBIS proj 5351

# TESTED for perl (v5.18.4)


# script to filter and normalise RSL/UMI counts in input libraries ONLY
# it is used to create a reference set of barcodes present in two tech replicates on the input libraries
# to be used for filtering "ghost" barcode-guide combinations in screen results

# RSL and UMI are used interchangeably

# based on earlier normaliseUMIcounts.input.v0.9.pl with the following changes


# substitute RSL for UMI in stdout messages
# better stdout messages
# use counts per million CPM rather than RPM


# input file should contain ONLY these columns
# column "Unknown" if present must be removed
# 

# RSL.guide,guide.set,Input1,Input2
# GM2A_4_TCTTAATAGA,GM2A_4,1,0
# GM2A_4_GGCGAGTTAT,GM2A_4,0,1
# GM2A_4_AACATGTTCC,GM2A_4,0,0
# GM2A_4_TTAACGTGGA,GM2A_4,4,3
# GM2A_4_CATATACAAC,GM2A_4,0,1

# RSL-guide filtering 
# read cutoff for RSL-guide combo taken from CLI parameter (mandatory)
# --CO all=N cutoff =N reads in *both* (all) replicates


# normalisation: counts per million CPM
# calculated AFTER the filtering based on number of raw reads
# only the normalisation factors are calculated, normalised reads are not output! this is for comparison of different filtering cutoffs

# calc CPM  norm factors
# norm factors will be cpm - counts per million
# n - reads in the sample
# n_tot = all reads in the library
# nf = 1e6/n_tot
# n_norm  = n * nf


# OUTPUT (main)
# added counts to each RSL/UMI-sgRNA combination which passed filtering based on set CO
# this is a master file used for processing screen results

# OUTPUT (diagnostic)
# output of frequencies of read counts to each RSL/UMI-sgRNA combination (raw counts after filtering) in each sample
# the counts in each tech replicate are added
# for histogram plotting

# OUTPUT (log)
# file with the following information for the report (tab delimited)
# cutoff - number of retained combinations (number of printed lines) - library size for each replicates (colSums)

## to do
# ignore sample "Unknown"



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use List::Util qw( sum );
use File::Basename;
use File::Path qw(make_path);
#use Data::Dumper qw(Dumper);
use List::SomeUtils qw( all any);
#use IO::Tee;
#use IO::File;


my $script_name="filter_RSL_input.v0.11.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/SAMPLE.UMIcounts.csv\n";
	print "--outdir: /path/to/outdir\n";
	print "--pref: prefix for output files\n";
	print "--CO: selected read cutoff for RSL-guide combination to be included\n";
	print "--CO all=N cutoff =N reads in *both* (all) replicates\n";
	print "-h prints this message\n";
}


else{
	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outdir=s'		=>	\(my $outdir),
		'pref=s'		=>	\(my $prefixout),
		'CO=s'		=>	\(my $FO_cli)
	) or die "Error in command line arguments";


	unless (-d $outdir) {
    	make_path($outdir);
	}


	#my @suffixlist=(".tsv",".csv",".txt");
	#my $basename = basename($infile,@suffixlist);

	my $basename=$prefixout;
	# all of these are output using the tech replicate cutoff $filt_CO_tech
	my $outfile_filt_fname="$basename\.filtered.csv"; #NOT normalised reads added for technical replicates after filtering using the tech replicate cutoff
	my $outfile_freqs_raw_fname="$basename\.frequencies\.raw_reads_aggregated.tsv"; #frequencies for aggregated reads after filtering using the tech replicate cutoff
	#my $outfile_norm_fname="$basename\.rpm.csv"; #rpm normalised reads added for technical replicates after filtering
	my $log_fname="$basename\.readme.log"; #readme with cmd and stdout echo
	my $log_report_fname="$basename\.report.log"; #readme with summary stats

	#my $outfile_norm="$outdir\/$outfile_norm_fname";
	my $outfile_filt="$outdir\/$outfile_filt_fname";
	my $outfile_freqs_raw="$outdir\/$outfile_freqs_raw_fname";
	my $log="$outdir\/$log_fname";
	my $log_report="$outdir\/$log_report_fname";


	open (LOG, ">", $log) or die "Cannot open logfile file $log: $!"; 

	open (STDOUT, "| tee $log") or die "Teeing off: $!\n";


	print "Project name is $basename\n\n";
	print "All output files $basename `.*` are of reads filtered based on filter $FO_cli\n\n";

	print "filtering parameters: $parameters\n\nfilter for technical replicates is $FO_cli (min number of reads in each of the replicates for each RSL-guide combination to be retained)\n\n";

	print "output files are \n$outfile_filt_fname - NOT normalised reads added for technical replicates after filtering using the cutoff $FO_cli\n";
	print "$outfile_freqs_raw_fname - frequencies for aggregated reads after filtering using the selected cutoff\n";
	print "$log_fname - log file which includes this text\n\n";

	open (LOG_REP, ">", $log_report) or die "Cannot open logfile file $log_report: $!"; 
	print LOG_REP "$FO_cli\t";



	#for calculation of normalisation factors
	my @col_sums;
	my $number_of_samples;

	#frequencies
	my %freq_aggreg_raw_reads; #hash where the keys are readcounts (aggregated) and the value is the frequency (for histogram plotting)

	#outfile header
	my $header_guide_barcode;
	my $header_smpls;
	my $header_per_guide_raw;


	#red the first line only to determine whether Unknown is present

	# open my $file, '<', "filename.txt"; 
	# my $firstLine = <$file>; 
	# close $file;


	print "Filter reads and calculate CPM normalisation factors\n";

	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 

	open (OUT_RAW_COUNTS_AGGREG, ">","$outfile_filt") or die "Cannot open output file $outfile_filt: $!";
	print OUT_RAW_COUNTS_AGGREG "RSL_guide,guide_set,InputFilt\n";

	#line counter
	my $line_cnt=0;

	while(<INFILE>){

		chomp $_;
		my @line=split /,/;


		unless ($_ =~m/RSL\.guide/){

			my $RSL_guide=shift @line;
			my $guide_set=shift @line;

			# if($_ =~m/Unknown/){ #?always the last column > assume this for now; remove this last element if Unknown has been matched
			# 	my $unknown_smpl=pop @line;
			# }
		
			$number_of_samples=scalar(@line);

			# apply filter for the presence of the read in all and any replicates before anything happens
  			if ( (all { $_ >= $FO_cli }(@line) ) && (any { $_ > 0 }(@line) ) ) {

  				#total read counts per sample $colsums
  				#@values_raw_readcounts is for 
				foreach my $i (0.. $#line) {
					$col_sums[$i] += $line[$i];
					#print "$i\t $col_sums[$i] \t $line[$i] \n";			
				}


				#aggregated count / for frequencies output
				my $aggreg_count=sum(@line);
				if (exists($freq_aggreg_raw_reads{$aggreg_count})){
					$freq_aggreg_raw_reads{$aggreg_count}++;
				}else{
					$freq_aggreg_raw_reads{$aggreg_count}=1;
				}

				my $line = "$RSL_guide,$guide_set,$aggreg_count";
				print OUT_RAW_COUNTS_AGGREG "$line\n";
				++$line_cnt;

			}

		# #get header
		}else{ #get header
			my $header1=shift @line;
			my $header2=shift @line;

			# if($_ =~m/Unknown/){ #?always the last column > assume this for now; remove this last element if Unknown has been matched
			# 	my $unknown_smpl=pop @line;
			# }

			my @UMI_counts=@line;

			$header_smpls=join ",", @UMI_counts;

			my $line_out="$header1,$header2,$header_smpls";

			$header_smpls=~s/,/\t/g;
			$header_per_guide_raw="$header1\tgene\t$header_smpls";

		}
	}
	print LOG_REP "$line_cnt\t";


	close(INFILE);
	close(OUT_RAW_COUNTS_AGGREG);
	
	print "Col sums are @col_sums\n";

	my $col_sums_log=join ",",@col_sums;
	print LOG_REP "$col_sums_log\n";
 

	my $tot_sum=sum(@col_sums);


	my $AvgColSum=$tot_sum/$number_of_samples;

	print "Total sum is $tot_sum\nAverage sum is $AvgColSum\tNumber of samples is $number_of_samples\n";

	my @CPM_nfs= map {1000000 / $_} @col_sums;
	my $cpm_nfs=join "\t",@CPM_nfs;

	my $aggreg_CPM_nf=1000000/$tot_sum;

	print "CPM normalisation factors are $cpm_nfs and aggregate CPM normalisation factor is $aggreg_CPM_nf\n\n";

	print "Print the output file\n";


	#aggreg frequencies
	open(OUTFILE_AGGREG_FREQ_RAW, ">","$outfile_freqs_raw") or die "Cannot open output file $outfile_freqs_raw: $!";

	for my $readcount (sort {$a <=> $b} keys %freq_aggreg_raw_reads){
		print OUTFILE_AGGREG_FREQ_RAW "$readcount\t$freq_aggreg_raw_reads{$readcount}\n";
	}

	close(OUTFILE_AGGREG_FREQ_RAW);

	#sanity checks
	#print "Array of values\n";
	#print "@values_raw_readcounts\n";

	#print "Dumper for the second hash\n";
	#print Dumper \%proc_freq_raw_reads;


	close(LOG_REP);
	close(LOG);
}
exit;

