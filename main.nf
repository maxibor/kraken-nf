#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
     kraken-nf: simple Kraken2/Bracken Nextflow pipeline
     Homepage: https://github.com/maxibor/kraken-nf
     Author: Maxime Borry <maxime_borry@eva.mpg.de>
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/kraken-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz' --krakendb '/path/to/minikraken2_v2_8GB_201904_UPDATE.tgz'
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --krakendb                    Path to MiniKraken2_v2_8GB Database
      --kraken_mem                  Memory allocated to kraken2. Ex: '4G'. Default to ${params.kraken_mem}

    Settings:
      --minimum_read_length         Minimum read length to keep. Default to ${params.minimum_read_length}
      --minhit                      Minimum number of Kraken hits to report Taxonomic level. Defaults to ${params.minhit}
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --pairedEnd                   Specified if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --build_bracken_db            Build Bracken database (true | false). Default = ${params.build_bracken_db}
      --run_bracken                 Run Bracken (true | false). Default = ${params.run_bracken}
      --bracken_level               Specifies the taxonomic level for Bracken. Default = ${params.bracken_level}
      --bracken_threshold           Specifies the threshold for Bracken. Default = ${params.bracken_threshold}

    Options:
      --results                     The output directory where the results will be saved. Defaults to ${params.results}
      --help  --h                   Shows this help page
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromFilePairs( params.reads, size: params.pairedEnd ? 2 : 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {reads_to_trim}

Channel
    .fromPath(params.krakendb, checkIfExists: true)
    .ifEmpty { exit 1, "Cannot find any KrakenDB matching: ${params.krakendb}\n" }
    .first()
    .set {krakendb}

process AdapterRemoval {
    tag "$name"

    label 'expresso'

    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq.gz') into trimmed_reads
        file("*.settings") into adapter_removal_results_multiqc

    script:
        settings = name+".settings"
        if (params.pairedEnd && !params.collapse){
            out1 = name+".pair1.trimmed.fastq.gz"
            out2 = name+".pair2.trimmed.fastq.gz"
            """
            AdapterRemoval \\
                --basename $name \\
                --file1 ${reads[0]} \\
                --file2 ${reads[1]} \\
                --trimns \\
                --trimqualities \\
                --minquality 20 \\
                --minlength ${params.minimum_read_length} \\
                --output1 $out1 \\
                --output2 $out2 \\
                --threads ${task.cpus} \\
                --qualitybase ${params.phred} \\
                --settings $settings
            """
        } else if (params.pairedEnd && params.collapse) {
            trimmed_out = name+".trimmed.fastq.gz"
            """
            AdapterRemoval  \\
                --qualitybase ${params.phred} \\
                --file1 ${reads[0]} \\
                --file2 ${reads[1]} \\
                --collapse \\
                --trimns \\
                --trimqualities \\
                --minquality 20 \\
                --minlength ${params.minimum_read_length} \\
                --basename $name \\
                --threads ${task.cpus} \\
                --settings $settings \\
                --seed 42 \\
                --gzip

                cat *.collapsed.gz *.collapsed.truncated.gz > $trimmed_out
            """
        } 
        else {
            se_out = name+".trimmed.fastq.gz"
            """
            AdapterRemoval \\
                --basename $name \\
                --file1 ${reads[0]} \\
                --minlength ${params.minimum_read_length} \\
                --trimns \\
                --trimqualities \\
                --minquality 20 \\
                --minlength ${params.minimum_read_length} \\
                --output1 $se_out \\
                --threads ${task.cpus} \\
                --qualitybase ${params.phred} \\
                --settings $settings
            """
        }       
}


process kraken2 {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/kraken", mode: 'copy', pattern: '*.kraken2_minimizer_report'

    input:
        set val(name), path(reads) from trimmed_reads
        path db from krakendb

    output:
        set val(name), file('*.kraken.out') into kraken_out
        set val(name), file('*.kraken2_minimizer_report') into kraken_report

    script:
        out = name+".kraken.out"
        kreport = name+".kraken2_minimizer_report"
        if (params.pairedEnd && !params.collapse){
            """
            kraken2 --db $db --threads ${task.cpus} --output $out --report-minimizer-data --report $kreport --paired ${reads[0]} ${reads[1]}
            """    
        } else {
            """
            kraken2 --db $db --threads ${task.cpus} --output $out --report-minimizer-data --report $kreport ${reads[0]}
            """
        }
        
}

kraken_report
    .into {kraken_report_parse; kraken_report_back} 

process kraken_report_backward_compatibility {
  tag "$prefix"

  label 'ristretto'

  publishDir "${params.results}/kraken", mode: 'copy', pattern: '*.kreport'

  input:
  tuple val(prefix), path(kraken_r) from kraken_report_back

  output:
  tuple prefix, path("*.kreport") into kraken_report_old
  script:
  kreport = prefix+".kreport"

  """
  cut -f1-3,6-8 $kraken_r > $kreport
  """
}

kraken_report_old
    .into {kraken_report_multiqc; kraken_report_bracken}

if (params.build_bracken_db) {
    process build_bracken_db {
    label 'expresso'

    publishDir "${params.results}/bracken_db", mode: 'copy', pattern: '*{kraken,kmer_distrib}'

    input:
        path db from krakendb
    output:
        path db_name into bracken_db
    script:
        db_name = db.toString()+"_bracken"
        """
        bracken-build \\
            -d $db \\
            -t $task.cpus \\
            -k ${params.kraken_kmer_len} \\
            -l ${params.minimum_read_length} \\

        mv $db $db_name
        cp $db_name/*.kraken . 
        cp $db_name/*.kmer_distrib .
        """
    }
} else {
    bracken_db = krakendb
}

if (params.run_bracken) {
    process bracken {
        tag "$prefix"

        label 'expresso'

        publishDir "${params.results}/bracken", mode: 'copy', pattern: '*bracken_report'

        input:
        tuple val(prefix), path(kraken_report) from kraken_report_bracken
        path db from bracken_db

        output:
        tuple prefix, path("*.bracken_report") into bracken_report
        tuple prefix, path("*.new_bracken_report") into bracken_new_report

        script:
        bracken_report = prefix+".bracken_report"
        bracken_new_report = prefix+".new_bracken_report"
        """
        bracken \\
            -d $db \\
            -i $kraken_report\\
            -o $bracken_report \\
            -w $bracken_new_report \\
            -r ${params.minimum_read_length} \\
            -l ${params.bracken_level} \\
            -t ${params.bracken_threshold} \\
        """
    }
}

process kraken_parse {
    tag "$name"

    label 'ristretto'

    errorStrategy 'ignore'

    publishDir "${params.results}/kraken", mode: 'copy'

    input:
        set val(name), file(kraken_r) from kraken_report_parse

    output:
        file('*_kraken_parsed.csv') into kraken_parsed

    script:
        read_out = name+".read_kraken_parsed.csv"
        kmer_out =  name+".kmer_kraken_parsed.csv"
        """
        kraken_parse.py -c ${params.minhit} -or $read_out -ok $kmer_out $kraken_r
        """
}

process kraken_merge {

    label 'ristretto'

    publishDir "${params.results}", mode: 'copy'

    input:
        file(csv_count) from kraken_parsed.collect()

    output:
        file('*.csv') into kraken_merged

    script:
        read_out = "kraken_read_count.csv"
        kmer_out = "kraken_kmer_duplication.csv"
        """
        merge_kraken_res.py -or $read_out -ok $kmer_out
        """    
}

kraken_report_multiqc
    .map {it -> it[1]}
    .set {kraken_report_multiqc_file} 

process multiqc {

    publishDir "${params.results}/multiqc", mode: 'copy'

    label 'ristretto'

    input:
        path('adapterRemoval/*') from adapter_removal_results_multiqc.collect().ifEmpty([])
        path('kraken/*') from kraken_report_multiqc_file.collect().ifEmpty([])
    output:
        path('*multiqc_report.html')
    script:
        """
        multiqc .
        """

}