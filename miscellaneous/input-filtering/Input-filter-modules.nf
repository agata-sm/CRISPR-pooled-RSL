// nextflow workflow for the crispr pipeline additional scripts for filtering input libraries
// modules for Input-filter.nf

// 8vi2022

// outdirs
params.filt="input_filtered"
params.filtOut="${params.outdir}/${params.filt}"

params.counterRSL="Counter"
params.counterRSLOut="${params.outdir}/${params.counterRSL}"

// scripts
params.scripts="${projectDir}/bin"


/// modules

process crispr_counter {
    publishDir params.counterRSLOut, mode:'copy'

    label 'big_mem'

    input:
    path "${params.projname}.properties"

    output:
    path "${params.projname}.csv", emit: countstable_ch
    path "counter.stdout.txt"
    path "counter.stdout.parsed.txt"


    script:
    """
    java -Xmx${task.memory.giga}g -jar ${params.crisprcounterpath}/CrisprCounter.jar ${params.properties} &> counter.stdout.txt

    perl ${params.scripts}/parseCrisprCounter.pl -i counter.stdout.txt -o counter.stdout.parsed.txt

    """

}


process filter_input {
    publishDir params.filtOut, mode:'copy'

    label 'small'

    input:
    tuple path(cntable), val(cutoff)

    output:
    path "${params.projname}.${cutoff}/*filtered.csv", emit: filt_input_list_ch
    path "${params.projname}.${cutoff}/*frequencies.raw_reads_aggregated.tsv"
    path "${params.projname}.${cutoff}/*readme.log"
    path "${params.projname}.${cutoff}/*report.log", emit: filt_input_log_ch

    path "${params.projname}.${cutoff}.${params.refdatapref}/*.RSL.perguide.tsv", emit: rsl_countstable_filt_ch
    path "${params.projname}.${cutoff}.${params.refdatapref}/*.frequencies.filt_reads.tsv"
    path "${params.projname}.${cutoff}.${params.refdatapref}/*.nonfilt_reads.perguide.tsv"
    path "${params.projname}.${cutoff}.${params.refdatapref}/*.readme.log"


    script:

    if ( "${params.usereference}"== "TRUE" ){

        """
        perl ${params.scripts}/filter_RSL_input.v0.12.pl --infile $cntable --pref ${params.projname} --outdir ${params.projname}.${cutoff} --CO $cutoff

        perl ${params.scripts}/processUMIcounts.v0.14.1.pl --filter CO=${params.filtRowSums} --pref ${params.refdatapref}.${cutoff} --infile ${params.refdatacnttable} --input_lib ${params.projname}.${cutoff}/${params.projname}.filtered.csv --outdir ${params.projname}.${cutoff}.${params.refdatapref} --input_lib_design $params.librarydesign

        """

    }else{

        """
        perl ${params.scripts}/filter_RSL_input.v0.12.pl --infile $cntable --pref ${params.projname} --outdir ${params.projname}.${cutoff} --CO $cutoff

        touch ${params.projname}.${cutoff}.${params.refdatapref}/NULL.RSL.perguide.tsv
        touch ${params.projname}.${cutoff}.${params.refdatapref}/NULL.frequencies.filt_reads.tsv
        touch ${params.projname}.${cutoff}.${params.refdatapref}/NULL.nonfilt_reads.perguide.tsv
        touch ${params.projname}.${cutoff}.${params.refdatapref}/NULL.readme.log

        """


    }


}



process report {
    publishDir params.outdir, mode:'copy'

    label 'small'

    input:
    path('*')

    output:
    path "report.${params.projname}"

   
    script:
    """
    cp -r ${params.projdir} .
    cp -r ${projectDir}/bin/report_template-input/* .
    Rscript input_report_launcher.R ${params.inputcnttable} ${params.filtOut} ${params.projname} ${params.refdatapref} ${params.usereference}
    """
}



