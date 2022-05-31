#!c:/perl/bin/perl.exe

#script to generate config file for CrisprCounter.jar

# input: path to fastqdir, samples.txt with file names for R1, path to library file to use, projectPrefix

# out: path to projectPrefix.properties for CrisprCounter.jar


use warnings;
use strict;
use diagnostics;
use Getopt::Long;

my $script_name="makeCounterConfig.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--template: /path/to/template.properties\n";
	print "--fastqdir: /path/to/dir/files_R1.fastq.gz\n";
	print "--samples: /path/to/samples.txt\n";
	print "--library: /path/to/library.csv e.g. Brunello_Library_USE_THIS_ONLY.csv\n";
	print "--prefix: PREF for file names\n";
	print "--outdir: /path/to/outdir where PREF.properties will be saved\n";

}else{

	my $parameters=join(' ', @ARGV);

	#commandline parsing for parameters
	GetOptions(
		'template=s'		=>	\(my $template),
		'fastqdir=s'		=>	\(my $fastqdir),
		'outdir=s'		=>	\(my $outdir),
		'samples=s'		=>	\(my $samples),
		'prefix=s'		=>	\(my $prefix),
		'library=s'		=>	\(my $library)
	) or die "Error in command line arguments";


	my $outfile="$outdir\/$prefix\.properties";
	#print "$outfile\n";
	my $out_cnttable="$prefix\.csv";

	# information on fastq files, R1
	my %fastq_files;
	my @samples_and_files;

	open (INFILE, "<","$samples") or die "Cannot open input file $samples: $!"; 

	while(<INFILE>){
		
		chomp $_;
		my @line=split /\t/;

		unless ($line[0] eq qw /file/){
			unless ($line[0] =~m/^\s*$/){
				if ($line[0]=~m/(\S+)_R1_001.fastq.gz/){
					my $fastqpref=$1;
					my $R1_read="$fastqdir\/$line[0]";
					my $R2_read="$fastqdir\/$fastqpref\_R2_001.fastq.gz";
					my $sample=$line[1];

					push @samples_and_files, [$sample,$R1_read,$R2_read];

				}else{print "entry $line[0] listed in $samples is not the expected format\n";}

			}
		}
	}
	close(INFILE);


	#final arrays
	my @sample_names;
	my @R1;
	my @R2;

	foreach my $i (0 .. $#samples_and_files){
			 	 		$sample_names[$i]=$samples_and_files[$i][0];
			 	 		$R1[$i]=$samples_and_files[$i][1];
			 	 		$R2[$i]=$samples_and_files[$i][2];
	}

	my $sample_names_string=join ",",@sample_names;
	my $R1_string=join ",",@R1;
	my $R2_string=join ",",@R2;

	#print "$sample_names_string\n$R1_string\n$R2_string\n";

	#fill in the template

	open (OUTFILE, ">", $outfile) or die "Cannot open output file $outfile: $!"; 

	open (INFILE, "<","$template") or die "Cannot open input file $template: $!"; 

	while(<INFILE>){

		chomp $_;
		
		if($_ =~m /^design =/){
			print OUTFILE "design = $library\n";
		}

		elsif($_ =~m /^out =/){
			print OUTFILE "out = $out_cnttable\n";
		}

		elsif($_ =~m /^fastq =/){
			print OUTFILE "fastq = $R1_string\n";
		}

		elsif($_ =~m /^umifastq =/){
			print OUTFILE "umifastq = $R2_string\n";
		}
		elsif($_ =~m /^samples =/){
			print OUTFILE "samples = $sample_names_string\n";
		}		
		else{
			print OUTFILE "$_\n";
		}

	}
	close (OUTFILE);

}
exit;