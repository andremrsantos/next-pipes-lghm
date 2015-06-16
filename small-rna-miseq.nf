#! /usr/bin/env nextflow

/**
 * Small RNA analysis pipeline for Illumina Sequencing plataforms (e.g. MiSeq).
 * This pipeline includes all processing steps to be executed from raw read
 * pre-processing to transcript counting.
 *
 * Authors
 * André M. Ribeiro-dos-Santos <andremrsantos@gmail.com>
 */

/**
 * Defining some parameters to identify input reads and processing arguments.
 **/
fastq_regex = /^(.+)\.(fastq|fq)(\.gz|)$/
running_stamp = new Date().format('yyyyMMdd')

params.input = "./raw/*.fastq.gz"
params.genome = "/data/resources/references/human/hg19/star_genome"
params.annot = ""
params.label = "transcript"
params.adapter = "ATCTCGTATGCCGTCTTCTGCTTG"
params.min_size = 18
params.quality_cutoff = 20
params.max_mismatch = 3
params.min_match = 16
params.base_dir = "./"

/**
 * Building Working Enviroment...
 */

Channel.from('alignment', 'read', 'report')
  .subscribe {
    f = file("${params.base_dir}/${it}")
    if (f.exists() && !f.isDirectory()) {
        error "ERROR: Unable to create '${f}'..."
    } else if (!f.exists()) {
        f.mkdir()
    }
  }

logger_file = file("${params.base_dir}/small-rna.${running_stamp}.log")
logger_file.setText('')

logger = Channel.create().subscribe {
    logger_file << "${it}\n"
    log.info(it)
}

align_copy = Channel.create().subscribe { 
  logger << "Copying file '${it.getFileName()}'..."
  it.copyTo("${params.base_dir}/alignment/") 
}
reads_copy = Channel.create().subscribe { 
  logger << "Copying file '${it.getFileName()}'..."
  it.copyTo("${params.base_dir}/read/")
}
reptr_copy = Channel.create().subscribe { 
  logger << "Copying file '${it.getFileName()}'..."
  it.copyTo("${params.base_dir}/report/") 
}

logger << "Small RNA Pipeline (from LGHM)"
logger << "------------------------------"
logger << "running_at   : ${new Date()}"
logger << "working dir  : ${file(params.base_dir)}"
logger << "genome       : ${params.genome}"
logger << "annotation   : ${params.annot}"
logger << "minimal size : ${params.min_size}"
logger << "qual cutoff  : ${params.quality_cutoff}"
logger << "max mismatch : ${params.max_mismatch}"
logger << "min mismatch : ${params.min_match}"
logger << "samples      : "

/**
 * Loading input file list
 */
raw_reads = Channel.create()

Channel.fromPath(params.input)
      .ifEmpty { error "ERROR: No read found matching: '${params.input_reads}'" }
      .map { tuple( sampleName(it), it ) }
      .tap( raw_reads )
      .subscribe {
        logger << "\t${it[0]} => ${it[1]}"
      }

process filter_reads {
    input:
    val adapter from params.adapter
    val min_size from params.min_size
    val qual_cutoff from params.quality_cutoff
    set val(name), file(path) from raw_reads

    output:
    set val(name), file("${name}.cutted.fastq") into cut_reads
    file "${name}.cutted.fastq" into reads_copy
    stdout into logger

    script:
    fpath = file("${params.base_dir}/read/${name}.cutted.fastq")
    if (fpath.exists())
      """
      echo "------------------------------"
      echo "Loading previous result '${fpath}' ..."
      cp ${fpath} ./
      echo "------------------------------"
      """
    else 
      """
      echo "------------------------------"
      echo "Cutting and Trimming Reads..."
      echo "Sample: ${name}"
      cutadapt -a ${adapter} \
        -m ${min_size} \
        -q ${qual_cutoff} \
        -o ${name}.cutted.fastq \
        ${path}
      echo "------------------------------"
      """
}

process alignment {
    input:
    set val(sample), file(sample_path) from cut_reads

    output:
    set val(sample), file("${sample}.align.bam") into alignments
    file "${sample}.align.bam" into align_copy
    stdout into logger

    script:
    """
    echo "------------------------------"
    echo "Align Reads to Reference..."
    echo "Sample: ${sample}"
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${params.genome} \
        --readFilesIn ${sample_path} \
        --outFileNamePrefix ${sample} \
        --outFilterMismatchNmax ${params.max_mismatch} \
        --outFilterMatchNmin ${params.min_match} \
        --outFilterMultimapScoreRange 0 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --alignIntronMax 1 \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outSAMprimaryFlag AllBestScore
    samtools view ${sample}Aligned.out.bam -x NH | \
        grep -Pv 'HI:i:([2-9]|[1-9][0-9]+)' | \
        samtools view -bS - > ${sample}.align.bam
    echo "------------------------------"
    """
}

process quantify {
    input:
    set val(sample), file(sample_path)

    output:
    set val(sample), file("${sample}.${label}count.txt")
    file "${sample}.count.txt" into reptr_copy
    stdout into logger
    stderr into logger

    script:
    """
    echo "------------------------------"
    htseq-count -a 0 -f bam ${sample_path} ${params.annot} \
        > ${sample}.${params.label}.count.txt
    echo "------------------------------"
    """
}


def sampleName(Path input) {
    final name = input.getFileName().toString()
    def matcher = (name =~ fastq_regex)
    if (matcher.matches()) {
        return matcher[0][1]
    }
    return name
}
