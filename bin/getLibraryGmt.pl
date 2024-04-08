#!c:/perl/bin/perl.exe

# author: Agata Smialowska

# NBIS proj 5351

# TESTED for perl (v5.18.4)

# script to generate library files necessary for RSL part of the CRISPR analysis pipeline CRISPR-pooled-RSL
# starting from library definition file in csv format

# input
# library confing file, w/o header
# CON249_1,GTAAACTTTGTCTGGAGTAT,CON249

# output:
# gmt file for GSEA
# text files for mageck - with control guides and genes



use warnings;
use strict;
use diagnostics;
use Getopt::Long;


my $script_name="getLibraryGmt.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/LibraryDefinition.csv in format sgRNA,sequence,gene\n";
	print "--outfile: /path/to/outfile.gmt\n";
	print "--outfile_con: /path/to/outfile.sgRNAcontrols.txt\n";	
	print "--outfile_gcon: /path/to/outfile.genecontrols.txt\n";	
	print "--outfile_lib: /path/to/outfile.librarydefinition.txt with corrected line endings (unix)\n";	
	print "-h prints this message\n";
}


else{
	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $input_lib_design),
		'outfile=s'		=>	\(my $outfile),
		'outfile_con=s'		=>	\(my $outfile_ctrl),
		'outfile_gcon=s'		=>	\(my $outfile_gctrl),
		'outfile_lib=s'		=>	\(my $outfile_lib)
	) or die "Error in command line arguments";

	open (OUTFILE_CTRL, ">", $outfile_ctrl) or die "Cannot open outfile $outfile_ctrl: $!";
	open (OUTFILE_GCTRL, ">", $outfile_gctrl) or die "Cannot open outfile $outfile_gctrl: $!";

	open(OUTFILE_LIB, ">", $outfile_lib) or die "Cannot open outfile $outfile_lib: $!";

	my %gene_guide;
	open (INFILE_INPUT_LIBDES, "<", $input_lib_design) or die "Cannot open file with processed input library $input_lib_design: $!";
	while (<INFILE_INPUT_LIBDES>){

    	# change line endings to \n
		# my $newline;
		# for ("\r", "\n", "\r\n", "\n\r", "\r\r", "\n\n") {
    	# 	$newline = 'text'.$_."text";
    	# 	$newline =~ s/$1/\n/g if m/(\r\n?|\n\r?)/;
    	# 	print OUTFILE_LIB "$newline\n";
    	# }
    	# chomp $newline;

    	# change line endings to \n from \r\n
 		$_=~s/[\r\n]+//;
		chomp $_;
		print OUTFILE_LIB "$_\n";

		my @line=split /,/;	

		unless ($_ =~m/Guide\t.*\tGene/){

			my $guide=$line[0];
			my $gene=$line[2];

			push(@{ $gene_guide{$gene} }, $guide);

			#$guide_gene{$guide}=$gene;

			if($guide =~m/CON\d+/){
				print OUTFILE_CTRL "$guide\n";
			}
		}

	}	
	close(INFILE_INPUT_LIBDES);
	close(OUTFILE_CTRL);

	close(OUTFILE_LIB);


	open (OUTFILE, ">", $outfile) or die "Cannot open outfile $outfile: $!";

	foreach my $gene (sort keys %gene_guide) {
		my @guides=@{ $gene_guide{$gene} };
		my $guides=join "\t",@guides;
        print OUTFILE "$gene\tna\t$guides\n";

        if($gene =~m/CON\d+/){
			print OUTFILE_GCTRL "$gene\n";
        }

    }

	close(OUTFILE);
	close(OUTFILE_GCTRL);

}

exit;

