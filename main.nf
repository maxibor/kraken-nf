#!/usr/bin/env nextflow




params.reads = ''
params.krakendb = '/path/to/minikraken2_v2_8GB_201904_UPDATE.tgz'
params.phred = 33
params.results = './results'
params.pairedEnd = true
params.help == false

def helpMessage() {
    log.info"""
     kraken-nf: simple Kraken2 Nextflow pipeline
     Homepage: https://github.com/maxibor/kraken-nf
     Author: Maxime Borry <borry@shh.mpg.de>
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/kraken-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz' --krakendb '/path/to/minikraken2_v2_8GB_201904_UPDATE.tgz'
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --krakendb                    Path to MiniKraken2_v2_8GB Database

    Settings:
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --pairedEnd                   Specified if reads are paired-end (true | false). Default = ${params.pairedEnd}

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


process AdapterRemoval {
    tag "$name"

    conda 'bioconda::adapterremoval'

    label 'expresso'

    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq') into trimmed_reads
        file("*.settings") into adapter_removal_results

    script:
        out1 = name+".pair1.trimmed.fastq"
        out2 = name+".pair2.trimmed.fastq"
        se_out = name+".trimmed.fastq"
        settings = name+".settings"
        if (params.pairedEnd){
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $out1 --output2 $out2 --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } else {
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $se_out --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        }
            
}

process miniKraken {
    tag "$name"

    conda 'bioconda::kraken2'

    label 'intenso'

    input:
        set val(name), file(reads) from trimmed_reads

    output:
        set val(name), file('*.kraken.out') into kraken_out
        set val(name), file('*.kreport') into kraken_report

    script:
        out = name+".kraken.out"
        kreport = name+".kreport"
        if (params.pairedEnd){
            """
            kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport --paired ${reads[0]} ${reads[1]}
            """    
        } else {
            """
            kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport ${reads[0]}
            """
        }
        
}

process kraken_parse {
    tag "$name"

    conda 'python=3.6'

    label 'ristretto'

    input:
        set val(name), file(kraken_r) from kraken_report

    output:
        set val(name), file('*.kraken_parsed.csv') into kraken_parsed

    script:
        out = name+".kraken_parsed.csv"
        """
        kraken_parse.py $kraken_r
        """    
}

process kraken_merge {

    conda 'python=3.6 pandas numpy'

    label 'ristretto'

    publishDir "${params.results}", mode: 'copy'

    input:
        file(csv_count) from kraken_parsed.collect()

    output:
        file('kraken_otu_table.csv') into kraken_merged

    script:
        out = "kraken_otu_table.csv"
        """
        merge_kraken_res.py -o $out
        """    
}