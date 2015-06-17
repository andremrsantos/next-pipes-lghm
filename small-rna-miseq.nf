#! /usr/bin/env nextflow

/**
 * Small RNA analysis pipeline for Illumina Sequencing plataforms (e.g. MiSeq).
 * This pipeline includes all processing steps to be executed from raw read
 * pre-processing to transcript counting.
 *
 * Authors
 * André M. Ribeiro-dos-Santos <andremrsantos@gmail.com>
 */

// Workflow parameters.
params.input          = "./raw/*.fastq.gz"
params.genome         = "/data/resources/references/human/hg19/star_genome"
params.annot          = "/data/resources/mirbase/20/hsa.gff"
params.feature        = "miRNA"
params.feature_id     = "Name"
params.label          = "transcript"
params.adapter        = "ATCTCGTATGCCGTCTTCTGCTTG"
params.min_size       = 18
params.quality_cutoff = 20
params.max_mismatch   = 3
params.min_match      = 16
params.base_dir       = "./"

// Setting workspace environment

fastq_regex   = /^(.+)\.(fastq|fq)(\.gz|)$/
running_stamp = new Date().format('yyyyMMdd')
read_folder   = "reads"
align_folder  = "alignments"
report_folder = "report"

logger_file   = file("${params.base_dir}/small-rna.${running_stamp}.log")
logger_file.setText('')

final logger = { String text ->
  logger_file << "$text\n"
  log.info(text)
}

final copy = { String folder, Path result ->
  destiny_folder = getFolder(folder)
  if( ! file("$destiny_folder/${result.getFileName()}").exists() )
    result.copyTo(destiny_folder)
}

// Building folder structure
Channel.from(read_folder, align_folder, report_folder)
  .subscribe {
    f = file( getFolder(it) )
    if (f.exists() && !f.isDirectory())
        error "ERROR: Unable to create '${f}'..."
    else if (!f.exists())
        f.mkdir()
  }

// Loading Input Path
raw_reads = Channel.create()
file_list = Channel.fromPath(params.input)
  .filter( ~fastq_regex )
  .ifEmpty { error "ERROR: No read found matching: '${params.input_reads}'" }
  .map { tuple( sampleName(it), it ) }
  .tap( raw_reads )

// Logging running informations
logger("Small RNA Pipeline (from LGHM)")
logger("------------------------------")
logger("running_at   : ${new Date()}")
logger("working dir  : ${file(params.base_dir)}")
logger("genome       : ${params.genome}")
logger("annotation   : ${params.annot}")
logger("minimal size : ${params.min_size}")
logger("qual cutoff  : ${params.quality_cutoff}")
logger("max mismatch : ${params.max_mismatch}")
logger("min mismatch : ${params.min_match}")
logger("feature      : ${params.feature}")
logger("feature_id   : ${params.feature_id}")
logger("samples      : ")
file_list.subscribe { name, path -> logger("  --> $name @ '$path'") }

// 1st STEP: QC
// Filtering and Trimming low quality reads and adpater contaminant removal
process filter_reads {
    input:
    val read_folder
    set val(name), file(path) from raw_reads

    output:
    set val(name), file("${name}.cutted.fastq") into cut_reads
    set val(read_folder), file("${name}.cutted.fastq") into qc_copy
    stdout into qc_logger

    script:
    rescueScript("${name}.cutted.fastq", read_folder, """
      echo "...."
      echo "Cutting and Trimming Reads..."
      echo "Sample: ${name}"
      cutadapt -a ${params.adapter} \
        -m ${params.min_size} \
        -q ${params.quality_cutoff} \
        -o ${name}.cutted.fastq ${path}
      """ )
}

qc_copy.subscribe(copy)
qc_logger.subscribe(logger)

// 2nd - Alignmento to Reference Genome
// The remaining reads are aligned to the reference genome using
// STAR aligner.
process alignment {
    input:
    val align_folder
    set val(sample), file(sample_path) from cut_reads

    output:
    set val(sample), file("${sample}.align.bam") into alignments
    set val(align_folder), file("${sample}.align.bam") into align_copy
    stdout into align_logger

    script:
    rescueScript( "${sample}.align.bam", align_folder, """
      echo "Aligning reads to reference ..."
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
      samtools view -h ${sample}Aligned.out.bam -x NH | \
        grep -Pv 'HI:i:([2-9]|[1-9][0-9]+)' | \
        samtools view -bS - > ${sample}.align.bam
      """ )
}

align_copy.subscribe(copy)
align_logger.subscribe(logger)

// 3rd - Quantify
// Transcripts are quantified using htseq-count
process quantify {
    input:
    val report_folder
    set val(sample), file(sample_path) from alignments

    output:
    set val(report_folder), file("${sample}.${params.label}.count.txt") into quant_copy
    stdout into quant_logger

    script:
    rescueScript("${sample}.${params.label}.count.txt", report_folder, """
      echo "....."
      echo "Quantifying transcripts ..."
      echo "Sample: ${sample}"
      htseq-count -a 0 -f bam \
        -t ${params.feature} -i ${params.feature_id} \
        ${sample_path} ${params.annot} \
        > ${sample}.${params.label}.count.txt
      """ )
}

quant_copy.subscribe(copy)
quant_logger.subscribe onNext: logger, onComplete: { logger("Concluded processing at: ${new Date()}") }

def sampleName(Path input) {
    final name = input.getFileName().toString()
    def matcher = (name =~ fastq_regex)
    if (matcher.matches()) {
        return matcher[0][1]
    }
    return name
}

def getFolder(String folder) {
  return "${params.base_dir}/$folder"
}

def rescueScript(name, folder, script) {
  final destiny_file = file("${getFolder(folder)}/$name")

  if ( destiny_file.exists() )
    """
    echo "Loading previous result @ '$destiny_file' ..."
    cp $destiny_file ./
    """
  else
    script
}
