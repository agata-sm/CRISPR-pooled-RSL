/* 
 * Pipeline for processing and reporting results of barcoded pooled CRISPR screens
 * Written for CRISPR Genomics Facility, SciLifeLab, Stockholm, Sweden
 * 
 * Author: Agata Smialowska
 * March - December 2022
 */ 


// modules for crispr-pooled-rsl.nf

// 17iii2022 - 19xii2022

// outdirs
params.readscnt="reads/count_table_mageck"
params.readsCntOut="${params.outdir}/${params.readscnt}"

params.readsrra="reads/rra_mageck"
params.readsRraOut="${params.outdir}/${params.readsrra}"

params.counterRSL="RSL/Counter"
params.counterRSLOut="${params.outdir}/${params.counterRSL}"

params.filterRSL="RSL/count_table_filtered"
params.filterRSLOut="${params.outdir}/${params.filterRSL}"

params.mageckRSL="RSL/rra_mageck"
params.mageckRSLOut="${params.outdir}/${params.mageckRSL}"

params.fastqc="FastQC"
params.fastqcOut="${params.outdir}/${params.fastqc}"

params.libdir="library"
params.libdirOut="${params.projdir}/${params.libdir}"

// assets
params.countertemplate="${projectDir}/assets/template.properties"


// scripts
params.scripts="${projectDir}/bin"
params.crisprcounterpath="/proj/sllstore2017103/nbis5351/For_NBIS" // path to executable CrisprCounter.jar

// versions
params.verfile="software.versions"

/// modules



process prep_library_files {
    //publishDir params.libdirOut, mode:'copy'

    label 'small'

    input:
    //path ${params.librarydesign}
    path lib_ch

    output:
    path "library.gmt" , emit: lib_gmt_ch
    path "library.ctrl_sgRNAs.txt" , emit: lib_ctrls_sgRNA_ch
    path "library.ctrl_genes.txt" , emit: lib_ctrls_gene_ch
    path "library.definition.txt" , emit: lib_definition_ch
    path "${params.verfile}"

    script:
    """
    #module load perl_modules/5.18.4

    perl ${params.scripts}/getLibraryGmt.pl --infile ${lib_ch} --outfile library.gmt --outfile_con library.ctrl_sgRNAs.txt --outfile_gcon library.ctrl_genes.txt --outfile_lib library.definition.txt

    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process ** prep_library_files **" >>${params.verfile}
    perl -v >>${params.verfile}
    echo "getLibraryGmt.pl" >>${params.verfile}
    """

}


process cp_library_files_reads {
    publishDir params.libdirOut, mode:'copy'

    label 'small'
    
    input:
    path lib_ch
    path ctrls_gene_ch

    output:
    path "reads/*"


    script:

        if ( "${params.mageckCountNorm}"== "control" ){
        
        """
        mkdir -p reads
        cp ${params.ctrl_file} reads
        cp ${lib_ch} reads
        """

        }else{

        """
        mkdir -p reads
        cp ${lib_ch} reads
        """

        }

}


process mageck_count_reads {
    publishDir params.readsCntOut, mode:'copy'

    label 'mid_mem'

    input:
    path fastqr1_ch
    val smpls_ch
    path ctrls_sgRNA_ch
    path ctrls_gene_ch

    output:
    path "${params.projname}.count.txt" , emit: count_table_reads_mageck_raw_ch
    path "${params.projname}.count_normalized.txt" , emit: count_table_reads_mageck_norm_ch
    path "${params.projname}.countsummary.txt" 
    path "${params.projname}.log"
    //path "${params.projname}*.R" // the "." in prefix are subbed with "_" in these files
    //path "${params.projname}*.Rnw"
    path "${params.projname}*.pdf"
    path "${params.verfile}"


    script:

    if ( "${params.mageckCountNorm}"== "control" ){

        """
        mageck count --norm-method ${params.mageckCountNorm} ${params.ctrl_type} ${params.ctrl_file} --pdf-report -l ${params.librarydesign} -n ${params.projname} --fastq ${fastqr1_ch} --sample-label ${smpls_ch}

        echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
        date >>${params.verfile}
        echo "process **  mageck_count_reads **" >>${params.verfile}
        echo "mageck" >>${params.verfile}
        mageck -v >>${params.verfile}

        """

    }else{

        """
        #module load bioinfo-tools
        #module load MAGeCK/0.5.9.4
        #module load R_packages/4.1.1
        #module load pandoc/2.17.1.1

        mageck count --norm-method ${params.mageckCountNorm} --pdf-report -l ${params.librarydesign} -n ${params.projname} --fastq ${fastqr1_ch} --sample-label ${smpls_ch}
        
        echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
        date >>${params.verfile}
        echo "process **  mageck_count_reads **" >>${params.verfile}
        echo "mageck" >>${params.verfile}
        mageck -v >>${params.verfile}
        """
    }

}

process mageck_rra_reads {
    publishDir params.readsRraOut, mode:'copy'

    label 'mid_mem'

    tag {comparisonID}

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
    path "${params.verfile}"


    script:
    """
    mkdir $comparisonID
    #module load bioinfo-tools
    #module load MAGeCK/0.5.9.4
    #module load R_packages/4.1.1
    #module load pandoc/2.17.1.1

    mageck test -k ${cnttable} -c ${smplRef} -t ${smplTreat} -n ${comparisonID}/${comparisonID} --norm-method none --pdf-report
    
    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process ** mageck_rra_reads **" >>${params.verfile}
    echo "mageck" >>${params.verfile}
    mageck -v >>${params.verfile}
    """

}

process report_reads {
    publishDir params.projdir, mode:'copy'

    label 'mid_mem'

    input:
    path('*')
    path sampleInfo_ch
    path comparisonsInfo_ch

    output:
    path "report.reads"
    path "results/reads/rra_annotation"
    path "${params.verfile}"


    script:
    """
    #module load  R_packages/4.1.1
    #module load pandoc/2.10.1

    cp -r ${params.projdir} .
    cp -r ${projectDir}/bin/report_template/* .
    
    mkdir ${params.projname}/metadata
    cp ${params.sampleinfo} ${params.projname}/metadata
    cp ${params.comparisons} ${params.projname}/metadata
  
      if [["${params.scatters}" != "none"]]
    do
        cp ${params.scatters} ${params.projname}/metadata
    done



    Rscript report_launcher.R ${params.projname} ${params.projname} reads ${params.organism} ${sampleInfo_ch} ${comparisonsInfo_ch} ${params.scatters}

    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process ** report_reads **" >>${params.verfile}
    R --version >>${params.verfile}
    echo "please check Session Info in the report for package versions" >>${params.verfile}
    """

}



process crispr_counter {
    publishDir params.counterRSLOut, mode:'copy'

    label 'big_mem'

    input:
    path fastqr1_ch

    output:
    path "${params.projname}.properties"
    path "${params.projname}.csv", emit: rsl_countstable_ch
    path "counter.stdout.txt"
    path "counter.stdout.parsed.txt"
    path "${params.verfile}"

    script:
    """
    #module load perl_modules/5.18.4

    perl ${params.scripts}/makeCounterConfig.pl --template $params.countertemplate --samples $params.sampleinfo --library $params.librarydesign --prefix $params.projname --outdir . --fastqdir $params.fastqdir
  
    java -Xmx${task.memory.giga}g -jar ${params.crisprcounterpath}/CrisprCounter.jar ${params.projname}.properties &> counter.stdout.txt

    perl ${params.scripts}/parseCrisprCounter.pl -i counter.stdout.txt -o counter.stdout.parsed.txt

    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process ** crispr_counter **" >>${params.verfile}
    perl -v >>${params.verfile}
    echo "makeCounterConfig.pl" >>${params.verfile}
    echo "parseCrisprCounter.pl" >>${params.verfile}
    echo "CrisprCounter.jar, no version specified for standard output" >>${params.verfile}
    java -jar /opt/myjar/CrisprCounter.jar >>${params.verfile}
    """

}

process cp_library_files_RSL {
    publishDir params.libdirOut, mode:'copy'

    label 'small'
    
    input:
    path lib_ch
    path lib_gmt_ch

    output:
    path "RSL/*"


    script:

    """
    mkdir -p RSL
    cp ${lib_gmt_ch} RSL
    cp ${lib_ch} RSL
    """
}


process filter_RSL {
    publishDir params.filterRSLOut, mode:'copy'

    label 'mid_mem'

    input:
    path rsl_countstable
    path lib_definition

    output:
    path "${params.projname}.RSL.perguide.tsv", emit: rsl_countstable_filt_ch
    path "${params.projname}.frequencies.filt_reads.tsv"
    path "${params.projname}.nonfilt_reads.perguide.tsv"
    path "${params.projname}.readme.log"
    path "${params.verfile}"

    script:
    """
    #module load perl_modules/5.26.2
        
    perl ${params.scripts}/processUMIcounts.v0.14.pl --filter CO=${params.filtRowSums} --infile ${rsl_countstable} --input_lib ${params.libraryinputfilt} --outdir . --input_lib_design ${lib_definition}

    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process **  filter_RSL **" >>${params.verfile}
    perl -v >>${params.verfile}
    echo "processUMIcounts.v0.14.pl" >>${params.verfile}
    """
    }



process mageck_rra_RSL {
    publishDir params.mageckRSLOut, mode:'copy'

    label 'mid_mem'

    tag {comparisonID}


    input:
    tuple path(cnttable), path(lib_gmt), val(comparisonID), val(smplRef), val(smplTreat)

    output:
    path "${comparisonID}/${comparisonID}.rank_log2FC.tsv"
    path "${comparisonID}/${comparisonID}.${params.projname}.log"
    path "${comparisonID}/${comparisonID}.${params.projname}.pathway_summary.txt"
    path "${comparisonID}/${params.projname}.RSL.perguide.log"
    path "${comparisonID}/${comparisonID}.${params.projname}.gene_rra_summary.txt", emit: rsl_rra_mageck_ch
    path "${comparisonID}/${comparisonID}.${params.projname}.readme"
    path "${params.verfile}"

    script:
    """
    #module load bioinfo-tools
    #module load MAGeCK/0.5.9.4
    #module load R_packages/4.1.1
    #module load pandoc/2.17.1.1
    #module load perl_modules/5.18.4

    mkdir $comparisonID
    perl ${params.scripts}/rank_log2FC.v0.3.pl -i ${cnttable} -o ${comparisonID}/${comparisonID}.rank_log2FC.tsv -r ${smplRef} -t ${smplTreat}

    mageck pathway --gmt-file ${lib_gmt} --method rra --ranking-column 4 --ranking-column-2 3 --gene-ranking ${comparisonID}/${comparisonID}.rank_log2FC.tsv -n ${comparisonID}/${comparisonID}.${params.projname}

    cp "${comparisonID}/${comparisonID}.${params.projname}.pathway_summary.txt" "${comparisonID}/${comparisonID}.${params.projname}.gene_rra_summary.txt"
    
    echo "file ${comparisonID}.${params.projname}.pathway_summary.txt is the original file output by mageck" >${comparisonID}/${comparisonID}.${params.projname}.readme
    echo "file ${comparisonID}.${params.projname}.gene_rra_summary.txt is identical to ${comparisonID}.${params.projname}.pathway_summary.txt and its name more precisely reflects its contents and provenance" >>${comparisonID}/${comparisonID}.${params.projname}.readme
   
    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process ** mageck_rra_RSL **" >>${params.verfile}
    echo "mageck" >>${params.verfile}
    mageck -v >>${params.verfile}
    perl -v >>${params.verfile}
    echo "rank_log2FC.v0.3.pl" >>${params.verfile}
    """
}



process report_RSL {
    publishDir params.projdir, mode:'copy'

    label 'mid_mem'

    input:
    path('*')
    path sampleInfo_ch
    path comparisonsInfo_ch

  
    output:
    path "report.RSL"
    path "results/RSL/rra_annotation"
    path "${params.verfile}"

    script:
    """
    #module load  R_packages/4.1.1
    #module load pandoc/2.10.1

    cp -r ${params.projdir} .
    cp -r ${projectDir}/bin/report_template/* .

    mkdir ${params.projname}/metadata
    cp ${params.sampleinfo} ${params.projname}/metadata
    cp ${params.comparisons} ${params.projname}/metadata

    Rscript report_launcher.R ${params.projname} ${params.projname} RSL ${params.organism} ${sampleInfo_ch} ${comparisonsInfo_ch} ${params.scatters}

    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process ** report_RSL **" >>${params.verfile}
    R --version >>${params.verfile}
    echo "please check Session Info in the report for package versions" >>${params.verfile}
    """

}

process fastqc {
    publishDir params.fastqcOut, mode:'copy'

    label 'small'

    input:
    path fastqr1

    output:
    path('*')
    path "${params.verfile}"

    script:
    """
    #module load bioinfo-tools
    #module load FastQC/0.11.9
    
    fastqc ${fastqr1}

    echo "Software versions for crispr-pooled-rsl.nf" >${params.verfile}
    date >>${params.verfile}
    echo "process ** fastqc **" >>${params.verfile}
    fastqc -v >>${params.verfile}
    """

}





