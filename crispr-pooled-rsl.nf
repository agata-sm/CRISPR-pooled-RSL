#! /usr/bin/env nextflow


/* 
 * 2nd take on crispr pipeline
 * Using Nextflow DSL2
 * 
 * Author: Agata Smialowska
 * March 2022
 */ 

nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.resdir = "results"
params.projdir = "$launchDir/${params.projname}"

params.outdir = "${params.projdir}/${params.resdir}"

params.logdir = 'logs'
params.metadatadir = 'metadata'

params.fastqR1="$params.fastqdir/*R1*fastq.gz"

// if input based filtering is not desired
params.libraryinputfilt=""


log.info """\
 CRISPR - N F   P I P E L I N E
 ===================================
 metadata: ${params.sampleinfo}
 fastq files directory: ${params.fastqdir}
 input library design : ${params.librarydesign}
 input library filtered : ${params.libraryinputfilt}

 outdir       : ${params.outdir}
 """
 .stripIndent()

println ""

/////////////////////////////
// process metadata files

// get the files and sample names
println "Samples"
filesf = file("$params.sampleinfo")
filesf.withReader {
    String line
    while( line = it.readLine() ) {
        println line

    }
}


// get the list of contrasts and samples for mageck
println "Comparisons"
comparisonf = file("$params.comparisons")
comparisonf.withReader {
    String line
    while( line = it.readLine() ) {
        println line

    }
}

println ""
println ""

/////////////////////////////
// channels

// samples channel
smpls_ch= Channel.fromPath(params.sampleinfo, checkIfExists:true)
	smpls_ch
	    .splitCsv(header:true, sep: '\t', strip: true)
	    .map{ (it.sample) }
	    .toList()
	    .toListString().replace(/[/,"").replace(/]/,"").replace(/ /,"")
	    //.view()
	    .set { smpls_ch }

// fastq file paths channel - list of paths
fastqr1_ch= Channel.fromPath(params.sampleinfo, checkIfExists:true)
	fastqr1_ch
	    .splitCsv(header:true, sep: '\t', strip: true)
	    .map{ (it.file) }
	    .collect { "${params.fastqdir}/$it" }
	    //.view()
	    .set { fastqr1_ch }


// fastq file paths channel - paths
fastqr1_ch2= Channel.fromPath(params.sampleinfo, checkIfExists:true)
	fastqr1_ch2
	    .splitCsv(header:true, sep: '\t', strip: true)
	    .map{ (it.file) }
	    .map { "${params.fastqdir}/$it" }
	    //.view()
	    .set { fastqr1_ch2 }


// comparisons
comparisons_ch= Channel.fromPath(params.comparisons, checkIfExists:true)
	comparisons_ch
		    .splitCsv(header:true, sep: '\t', strip: true)
		    .map{ row-> tuple(row.name, row.reference, row.treatment) }
		    //.view()
		    .set { comparisons_ch }


// library definition
lib_ch= Channel.fromPath(params.librarydesign, checkIfExists:true)
	lib_ch
		.view()
		.set { lib_ch }



/////////////////////////////
// processes
include { prep_library_files; mageck_count_reads; mageck_rra_reads; report_reads; crispr_counter; filter_RSL; mageck_rra_RSL; report_RSL; fastqc } from './crisprRSL-modules.nf'



/////////////////////////////
// workflows

//default reads

workflow {

	//prep library files
	prep_library_files(lib_ch)

	//count reads
	mageck_count_reads(fastqr1_ch, smpls_ch)

	// mageck contrasts RRA reads
	cntReads_ch=mageck_count_reads.out.count_table_reads_mageck_norm_ch
		cntReads_ch
			.combine(comparisons_ch)
			//.view()
			.set { cntReads_ch }

	mageck_rra_reads(cntReads_ch)

	//report
	mageck_res_reads_gene_ch=mageck_rra_reads.out.gene_summary_reads_ch
	report_reads(mageck_res_reads_gene_ch.collect())

	//QC
	fastqc(fastqr1_ch2)

}

//alternative RSL

workflow RSL {


//workflow {

	//prep library files
	prep_library_files(lib_ch)

	// count reads
	crispr_counter(fastqr1_ch)

	filter_RSL(crispr_counter.out.rsl_countstable_ch)

	// mageck contrasts RSL
	cntRSL_ch=filter_RSL.out.rsl_countstable_filt_ch
	 	cntRSL_ch
	 		.combine(prep_library_files.out.lib_gmt_ch)
	 		.combine(comparisons_ch)
	 		.view()
	 		.set { cntRSL_ch }

	mageck_rra_RSL(cntRSL_ch)

	// //report
	mageck_res_RSL_gene_ch=mageck_rra_RSL.out.rsl_rra_mageck_ch
	report_RSL(mageck_res_RSL_gene_ch.collect())

	//QC
	//fastqc(fastqr1_ch2)

}
