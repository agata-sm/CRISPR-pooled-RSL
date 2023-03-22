	FROM perl:5.18.4	
	
	LABEL crisprcounter-perl based on perl:5.18.4
	LABEL use for running input filtering pipeline, as part of the CRISPR-pooled-RSL pipeline

	WORKDIR .

	RUN cpanm install Getopt::Long
	RUN cpanm install List::Util
	RUN cpanm install File::Basename
	RUN cpanm install File::Path
	RUN cpanm install List::SomeUtils
	RUN cpanm install experimental


	WORKDIR ./app


