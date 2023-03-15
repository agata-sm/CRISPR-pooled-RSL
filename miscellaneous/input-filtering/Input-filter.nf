#! /usr/bin/env nextflow


/* 
 * Filtering the artifacts of input libraries sequenced in tech replicates
 * In preparation for use in the barcoded CRISPR screen analysis
 * Using Nextflow DSL2
 * 
 * Author: Agata Smialowska
 * May 2022
 */ 

nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.resdir = "results"
params.projdir = "$launchDir/${params.projname}"

params.outdir = "${params.projdir}/${params.resdir}"


log.info """\
 INPUT FILTERING - N F   P I P E L I N E
 =========================================
 input library count table : ${params.inputcnttable}

 outdir       : ${params.outdir}
 """
 .stripIndent()

println ""



/////////////////////////////
// channels

//cutoffs_ch= Channel.of( 3, 5, 7 )

cutoffs_ch= Channel.of( 0, 1, 2, 3, 4, 5, 7, 10 )

// cnts_ch= Channel.fromPath(params.inputcnttable, checkIfExists:true)
// 	cnts_ch
// 		.combine(cutoffs_ch)
// 		.view()
// 		.set {cnts_ch}


input_properties_ch=Channel.fromPath(params.properties, checkIfExists:true)


/////////////////////////////
// processes
include { crispr_counter; filter_input; report } from './Input-filter-modules.nf'


/////////////////////////////
// workflows

//default

workflow {

	crispr_counter(input_properties_ch)

	//filter input library using different cutoffs
	filter_input(crispr_counter.out.countstable_ch)

	//report
	input_filt_ch=filter_input.out.rsl_countstable_filt_ch
	
	// input_filt_ch_col=input_filt_ch.collect()
	// 	input_filt_ch_col
	// 		.view()

	report(input_filt_ch.collect())

}
