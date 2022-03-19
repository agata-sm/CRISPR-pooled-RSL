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


/////////////////////////////
// process metadata files

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

// get the files and sample names
println "Samples"
filesf = file("$params.sampleinfo")
filesf.withReader {
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
	   	////.collect()
	   	////.toString().replace(/[/,"").replace(/]/,"").replace(/ /,"")
	   	////.toString()
	    .toListString().replace(/[/,"").replace(/]/,"").replace(/ /,"")
	    //.view()
	    .set { smpls_ch }

// fastq file paths channel
fastqr1_ch= Channel.fromPath(params.sampleinfo, checkIfExists:true)
	fastqr1_ch
	    .splitCsv(header:true, sep: '\t', strip: true)
	    .map{ (it.file) }
	    ////.map{ "${params.fastqdir}"{"$it/*"}  }

	    ////.map{ (it.file, it.sample) }
	    .collect { "${params.fastqdir}/$it" }
	    ////.toList()
	    //.view()
	    .set { fastqr1_ch }


// comparisons
comparisons_ch= Channel.fromPath(params.comparisons, checkIfExists:true)
	comparisons_ch
		    .splitCsv(header:true, sep: '\t', strip: true)
		    .map{ row-> tuple(row.name, row.reference, row.treatment) }
		    //.view()
		    .set { comparisons_ch }




/////////////////////////////
// chunks
include { mageck_count_reads; mageck_rra_reads; report_reads } from './crisprRSL-modules.nf'


/////////////////////////////
// workflows

//default reads

workflow {

	// count reads
	mageck_count_reads(fastqr1_ch, smpls_ch)

	// mageck RRA reads
	cntReads_ch=mageck_count_reads.out.count_table_reads_mageck_norm_ch
		cntReads_ch
			.combine(comparisons_ch)
			//.view()
			.set { cntReads_ch }

	mageck_rra_reads(cntReads_ch)

	//report
	mageck_res_reads_gene_ch=mageck_rra_reads.out.gene_summary_reads_ch
	report_reads(mageck_res_reads_gene_ch.collect())


}

//alternative RSL



