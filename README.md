[![Build Status](https://travis-ci.org/maxibor/kraken-nf.svg?branch=master)](https://travis-ci.org/maxibor/kraken-nf)

# Kraken-nf

Simple Kraken2 Nextflow pipeline

## Dependancies

- [conda](https://conda.io/en/latest/)
- [Nextflow](https://www.nextflow.io/) : `conda install -c bioconda nextflow`

## Usage

```
nextflow run maxibor/kraken-nf --reads "/path/to/paired_end_reads_*.{1,2}.fastq.gz" --krakendb "/path/to/minikraken2_v2_8GB_201904_UPDATE.tgz"
```

### Input

#### --reads

Use this to specify the location of your input FastQ files. For example:

`--reads 'path/to/data/sample_*_{1,2}.fastq'`

**Please note the following requirements:**

- The path must be enclosed in quotes
- The path must have at least one \* wildcard character
- When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs.

#### --krakendb

Path to Kraken2 MiniKraken2_v2_8GB Database. Can be downloaded [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)

#### --kraken_mem

Amount of memory allocated to Kraken2.
Depends on the kraken database.

Example:

```bash
--kraken_mem '9G'
```

### Output

#### kraken_otu_table.csv

OTU table of all input samples.  
Samples in columns, TAXID in rows

## Kraken2 database

Any Kraken2 database can be used, but the _minikraken2_v2_8GB_ is a good compromise between speed and accuracy.  
Please have a look at the [Kraken2 downloads page](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads) to download it.

## Help

```
$ nextflow run maxibor/kraken-nf --help
N E X T F L O W  ~  version 19.04.0
Launching `maxibor/kraken-nf` [compassionate_morse] - revision: f534a6a703
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
  --kraken_mem                  Memory allocated to kraken2. Ex: '4G'. Default to 9G

Settings:
  --minhit                      Minimum number of Kraken hits to report Taxonomic level. Defaults to 50
  --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 64
  --pairedEnd                   Specified if reads are paired-end (true | false). Default = true

Options:
  --results                     The output directory where the results will be saved. Defaults to ./results
  --help  --h                   Shows this help page
```
