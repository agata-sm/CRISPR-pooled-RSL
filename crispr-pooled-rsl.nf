#! /usr/bin/env nextflow


/* 
 * Pipeline for processing and reporting results of barcoded pooled CRISPR screens
 * Written for CRISPR Genomics Facility, SciLifeLab, Stockholm, Sweden
 * 
 * Author: Agata Smialowska
 * March - December 2022
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

params.fastqR1 = "$params.fastqdir/*R1*fastq.gz" //// OBS! this has probably caused unexplained behaviour when more than one instance of R1 is present in file names; fix this!



// if input based filtering is not desired
params.libraryinputfilt=""

// library control files
if( "${params.mageckCountNorm}" == "control" ){
	
	if( "${params.mageckCountCtrl}" == "sgRNA"){
		params.ctrl_type="--control-sgrna"
		
		
		//if (typeof "${params.control_file}" !== 'undefined') {
		//if( "${params.control_file}" != null){
		if( "${params.control_file}" ){

			params.ctrl_file="${params.control_file}"
		}else{
			params.ctrl_file="library.ctrl_sgRNAs.txt"
		}
	}
	else if( "${params.mageckCountCtrl}" == "gene"){
		params.ctrl_type="--control-gene"

		//if (typeof "${params.control_file}" !== 'undefined') {
		if( "${params.control_file}" ){
			params.ctrl_file="${params.control_file}"
		}else{
			params.ctrl_file="library.ctrl_genes.txt"
		}

	}

	//if (typeof "${params.control_file}" !== 'undefined') {
	if( "${params.control_file}" ){
		params.libctrl_string="${params.ctrl_file}"
	}else{
		params.libctrl_string="${params.ctrl_file} containing features CON* from ${params.librarydesign}"
	}


}else{
	params.libctrl_string="n.a."
	params.ctrl_type="n.a."

}

// check if scatters.txt is defined



log.info """\
 CRISPR - N F   P I P E L I N E
 ===================================
 metadata: ${params.sampleinfo}
 fastq files directory: ${params.fastqdir}
 input library design : ${params.librarydesign}
 input library filtered : ${params.libraryinputfilt}

 read normalisation: ${params.mageckCountNorm}
 control features type: ${params.ctrl_type}
 control file: ${params.libctrl_string}

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
		//.view()
		.set { lib_ch }


// input library filtering
inputfilt_ch= Channel.fromPath(params.libraryinputfilt, checkIfExists:true)
	inputfilt_ch
		//.view()
		.set { inputfilt_ch }


//metadata channels
sampleInfo_ch=Channel.fromPath(params.sampleinfo, checkIfExists:true)
comparisonsInfo_ch=Channel.fromPath(params.comparisons, checkIfExists:true)


/////////////////////////////
// processes
include { prep_library_files; cp_library_files_reads;cp_library_files_RSL; mageck_count_reads; mageck_rra_reads; report_reads; crispr_counter; filter_RSL; mageck_rra_RSL; report_RSL; fastqc } from './crisprRSL-modules.nf'



/////////////////////////////
// workflows

//default reads

workflow {

	//prep library files
	prep_library_files(lib_ch)
	ctrls_sgRNA_ch=prep_library_files.out.lib_ctrls_sgRNA_ch
	ctrls_gene_ch=prep_library_files.out.lib_ctrls_gene_ch
	cp_library_files_reads(lib_ch, prep_library_files.out.lib_ctrls_gene_ch)

	//count reads
	mageck_count_reads(fastqr1_ch, smpls_ch, ctrls_sgRNA_ch, ctrls_gene_ch)

	// mageck contrasts RRA reads
	cntReads_ch=mageck_count_reads.out.count_table_reads_mageck_norm_ch
		cntReads_ch
			.combine(comparisons_ch)
			//.view()
			.set { cntReads_ch }

	mageck_rra_reads(cntReads_ch)

	//report
	mageck_res_reads_gene_ch=mageck_rra_reads.out.gene_summary_reads_ch
	report_reads(mageck_res_reads_gene_ch.collect(), sampleInfo_ch, comparisonsInfo_ch)

	//QC
	fastqc(fastqr1_ch2)

}

//alternative RSL

workflow RSL {


//workflow {
	
	//prep library files
	prep_library_files(lib_ch)
	cp_library_files_RSL(lib_ch, prep_library_files.out.lib_gmt_ch)

	// count reads
	crispr_counter(fastqr1_ch)

	filter_RSL(crispr_counter.out.rsl_countstable_ch, prep_library_files.out.lib_definition_ch, inputfilt_ch)

	// mageck contrasts RSL
	cntRSL_ch=filter_RSL.out.rsl_countstable_filt_ch
	 	cntRSL_ch
	 		.combine(prep_library_files.out.lib_gmt_ch)
	 		.combine(comparisons_ch)
	 		//.view()
	 		.set { cntRSL_ch }

	mageck_rra_RSL(cntRSL_ch)

	// //report
	mageck_res_RSL_gene_ch=mageck_rra_RSL.out.rsl_rra_mageck_ch
	report_RSL(mageck_res_RSL_gene_ch.collect(), sampleInfo_ch, comparisonsInfo_ch)

	//QC
	//fastqc(fastqr1_ch2)

}
