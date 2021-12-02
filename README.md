[![Build Status](https://travis-ci.org/maxibor/kraken-nf.svg?branch=master)](https://travis-ci.org/maxibor/kraken-nf)

# Kraken-nf

Simple [Kraken2](https://github.com/DerrickWood/kraken2)/[Bracken](https://github.com/jenniferlu717/Bracken) Nextflow pipeline

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

Path to Kraken2 MiniKraken2_v2_8GB Database. A pre-build database can be downloaded [here](https://benlangmead.github.io/aws-indexes/k2)

#### --kraken_mem

Amount of memory allocated to Kraken2.
Depends on the kraken database.

Example:

```bash
--kraken_mem '9G'
```

## Kraken2 database

Any Kraken2 database can be used, but the _PlusPFP-8_ is a good compromise between speed and accuracy.  
Please have a look at the [Index zone](https://benlangmead.github.io/aws-indexes/k2) to download it.

## Help

```
$ nextflow run maxibor/kraken-nf --help
N E X T F L O W  ~  version 21.04.0

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
  --kraken_mem                  Memory allocated to kraken2. Ex: '4G'. Default to 9G

Settings:
  --minimum_read_length         Minimum read length to keep. Default to null
  --minhit                      Minimum number of Kraken hits to report Taxonomic level. Defaults to 10
  --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 33
  --pairedEnd                   Specified if reads are paired-end (true | false). Default = true
  --build_bracken_db            Build Bracken database (true | false). Default = false
  --run_bracken                 Run Bracken (true | false). Default = true
  --bracken_level               Specifies the taxonomic level for Bracken. Default = S
  --bracken_threshold           Specifies the threshold for Bracken. Default = 10

Options:
  --results                     The output directory where the results will be saved. Defaults to ./results
  --help  --h                   Shows this help page
```
