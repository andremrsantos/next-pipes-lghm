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


params.input_reads = "$baseDir/raw/*.fastq.gz"
params.genome = "/data/resources/references/human/hg19/star_genome"
params.annot = ""
params.label = ""
params.adapter = "ATCTCGTATGCCGTCTTCTGCTTG"
params.min_size = 18
params.quality_cutoff = 20
params.max_mismatch = 3
params.min_match = 16

/**
 * Building Working Enviroment...
 */
Channel.from('alignment', 'read', 'report')
  .subscribe {
    f = file("${baseDir}/${it}")
    if (f.exists() && !f.isDirectory()) {
        error "Unable to create '${f}'..."
    } else if (!f.exists()) {
        f.mkdir()
    }
  }

logger_file = file("${baseDir}/small-rna.${running_stamp}.log")
logger = Channel.subscribe {
    logger_file << "${it}\n"
    log.info(it)
}

align_copy = Channel.subscribe { it.copyTo("${baseDir}/alignment/") }
reads_copy = Channel.subscribe { it.copyTo("${baseDir}/reads/") }
reptr_copy = Channel.subscribe { it.copyTo("${baseDir}/report/") }

logger << "Small RNA Pipeline (from LGHM)"
logger << "------------------------------"
logger << "running_at   : ${new Date()}"
logger << "dir          : ${baseDir}"
logger << "genome       : ${params.genome.baseName()}"
logger << "annotation   : ${params.annot.baseName()}"
logger << "minimal size : ${params.min_size}"
logger << "qual cutoff  : ${params.quality_cutoff}"
logger << "max mismatch : ${params.max_mismatch}"
logger << "min mismatch : ${params.min_match}"
logger << "samples      : "

/**
 * Loading input file list
 */
Channel.fromPath(params.input_reads)
        .filter { it.FileName() =~ fastq_regex }
        .ifEmpty { error "No read found matching: '${params.input_reads}'" }
        .map { path -> tuble(readPrefix(path), path) }
        .each {
    logger << "\t${name} => ${path}"
}
.into(raw_reads)

if (param.label != '')
    label = "${params.label}-"
else
    label = ''

process filter_reads {
    input:
    set val(sample), file(sample_path) from raw_reads

    output:
    set val(sample), file("${sample}.cutted.fastq") into cut_reads
    file "${sample}.cutted.fastq" into reads_copy
    stdout into logger
    stderr into logger

    script:
    """
    echo "------------------------------"
    echo "Cutting and Trimming Reads..."
    echo "Sample: ${sample}"
    cutadapt -a ${params.adapter} \
        -m ${params.min_size} \
        -q ${params.quality_cutoff} \
        -o ${sample}.cutted.fastq \
        ${sample_path}
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
    stderr into logger

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
        samtools view -bS > ${sample}.align.bam
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
        > ${sample}.${label}count.txt
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
