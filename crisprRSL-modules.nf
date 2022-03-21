// second try of a nextflow workflow for the crispr pipeline
// modules for crispr-pooled-rsl.nf

// 17iii2022

// outdirs
params.readscnt="reads/count_table_mageck"
params.readsCntOut="${params.outdir}/${params.readscnt}"

params.readsrra="reads/rra_mageck"
params.readsRraOut="${params.outdir}/${params.readsrra}"

params.counterRSL="RSL/Counter"
params.counterRSLOut="${params.outdir}/${params.counterRSL}"

params.filterRSL="RSL/count_table_filtered"
params.filterRSLOut="${params.outdir}/${params.filterRSL}"


// assets
params.countertemplate="${projectDir}/assets/template.properties"


// scripts
params.scripts="${projectDir}/bin"

/// modules

process mageck_count_reads {
    publishDir params.readsCntOut, mode:'copy'

    input:
    path fastqr1_ch
    val smpls_ch

    output:
    path "${params.projname}.count.txt" , emit: count_table_reads_mageck_raw_ch
    path "${params.projname}.count_normalized.txt" , emit: count_table_reads_mageck_norm_ch
    path "${params.projname}.countsummary.txt" 
    path "${params.projname}.log"
    //path "${params.projname}_countsummary.R" // the "." in prefix are subbed with "_" in these files
    //path "${params.projname}_countsummary.Rnw"
    //path "${params.projname}_countsummary.pdf"

    script:
    """
    echo $smpls_ch

    mageck count --norm-method total --pdf-report -l $params.librarydesign -n $params.projname --fastq $fastqr1_ch --sample-label $smpls_ch
    """

}


process mageck_rra_reads {
    publishDir params.readsRraOut, mode:'copy'

    input:
    tuple path(cnttable), val(comparisonID), val(smplRef), val(smplTreat)

    output:
    path "${comparisonID}/${comparisonID}.sgrna_summary.txt"
    path "${comparisonID}/${comparisonID}.gene_summary.txt", emit: gene_summary_reads_ch
    path "${comparisonID}/${comparisonID}.R"
    path "${comparisonID}/${comparisonID}.log"
    path "${comparisonID}/${comparisonID}.pdf"
    path "${comparisonID}/${comparisonID}_summary.log"
    path "${comparisonID}/${comparisonID}_summary.Rnw"
    path "${comparisonID}/${comparisonID}_summary.pdf"


    script:
    """
    mkdir $comparisonID
    mageck test -k $cnttable -c $smplRef -t $smplTreat -n ${comparisonID}/${comparisonID} --norm-method none --pdf-report
    """

}

process report_reads {
    publishDir params.projdir, mode:'copy'

    input:
    path('*')

    output:
    path "report.reads"
    path "results/reads/rra_annotation"

    script:
    """
    #module load  R_packages/4.1.1
    #module load pandoc/2.10.1

    cp -r ${params.projdir} .
    mkdir ${params.projname}/metadata
    cp ${params.sampleinfo} ${params.projname}/metadata
    cp ${params.comparisons} ${params.projname}/metadata
    cp -r ${projectDir}/bin/report_template/* .
    Rscript report_launcher.R $params.projname $params.projname reads $params.organism
    """

}


process crispr_counter {
    publishDir params.counterRSLOut, mode:'copy'

    input:
    path fastqr1_ch

    output:
    path "${params.projname}.properties"
    path "${params.projname}.csv", emit: rsl_countstable_ch
    path "counter.stdout.summary.txt"

    script:
    """
    echo "$fastqr1_ch"
    #mkdir -p ${params.counterRSL}
    perl ${params.scripts}/makeCounterConfig.pl --template $params.countertemplate --samples $params.sampleinfo --library $params.librarydesignRSL --prefix $params.projname --outdir . --fastqdir $params.fastqdir
  
    # Rackham
    # correct the path
    #java -Xmx48G -jar /proj/snic2020-16-213/For_NBIS/CrisprCounter.jar ${params.counterRSL}/${$params.projname}.properties &> counter.stdout.txt

    #LOCAL tsts
    cp /Users/agata.smialowska/NBISproj/5351_CRISPR_pipeline/data/heldin_counter_stdout/counter.stdout.txt .
    cp ${params.cnttable} .

    perl ${params.scripts}/parseCrisprCounter.pl -i counter.stdout.txt -o counter.stdout.summary.txt

    """

}

process filter_RSL {
    publishDir params.filterRSLOut, mode:'copy'



}






