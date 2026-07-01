// ITS extraction for the fungal (marker==ITS) read-prep branch.
//
// 16S goes cutadapt -> DADA2 directly. ITS reads carry conserved rRNA flanks
// (18S / 5.8S / 28S) around the variable ITS1-5.8S-ITS2 region, so before DADA2
// we pull out the full ITS span with itsxrust (a fast Rust reimplementation of
// ITSx). Running on the cutadapt-trimmed reads keeps the same per-sample
// (sampleID, FASTQ) tuple flowing downstream — DADA2 then sees ITS-only reads,
// which is why the ITS length window (params.min_len/max_len) is short.
//
// itsxrust needs the `extract` subcommand and an explicit anchor HMM profile
// (--hmm; the bundled fungal F.hmm at params.itsx_hmm, staged into the task).
// `--preset hifi` selects the PacBio CCS-tuned HMM thresholds; `--derep` collapses
// identical extracted sequences; `--hmmer-cpu` sets the nhmmer threads; the full-ITS
// region is requested so the ITS1-5.8S-ITS2 span (not just ITS1 or ITS2) reaches DADA2.

process itsx_extract {
    conda (params.enable_conda ? "${projectDir}/envs/itsxrust.yml" : null)
    container params.itsx_container
    publishDir "${params.outdir}/itsx_extracted_FASTQ", pattern: '*.its.fastq.gz', mode: params.publish_dir_mode
    cpus { params.itsx_threads as int }
    label 'itsx'

    input:
    tuple val(sampleID), path(sampleFASTQ)
    path hmm

    output:
    tuple val(sampleID), path("${sampleID}.its.fastq.gz"), emit: fastq

    script:
    """
    itsxrust extract \\
        --input ${sampleFASTQ} \\
        --output ${sampleID}.its.fastq.gz \\
        --hmm ${hmm} \\
        --region full \\
        --preset hifi \\
        --derep \\
        --hmmer-cpu ${task.cpus}
    """
}
