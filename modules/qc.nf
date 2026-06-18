process QC_fastq {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"
    label 'cpu8'
    publishDir "$params.outdir/filtered_input_FASTQ", pattern: '*filterQ*.fastq.gz', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    path "${sampleID}.seqkit.readstats.tsv", emit: all_seqkit_stats
    path "${sampleID}.seqkit.summarystats.tsv", emit: all_seqkit_summary
    tuple val(sampleID), path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq
    path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq_files

    script:
    if (params.downsample > 0)
    """
    seqkit fx2tab -j $task.cpus -q --gc -l -H -n -i $sampleFASTQ |\
        csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.readstats.tsv
    seqkit stats -T -j $task.cpus -a ${sampleFASTQ} |\
        csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.summarystats.tsv
    seqkit seq -j $task.cpus --min-qual $params.filterQ $sampleFASTQ |\
        seqkit head -n $params.downsample --out-file ${sampleID}.filterQ${params.filterQ}.fastq.gz
    """
    else
    """
    seqkit fx2tab -j $task.cpus -q --gc -l -H -n -i $sampleFASTQ |\
        csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.readstats.tsv
    seqkit stats -T -j $task.cpus -a ${sampleFASTQ} |\
        csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.summarystats.tsv
    seqkit seq -j $task.cpus --min-qual $params.filterQ $sampleFASTQ --out-file ${sampleID}.filterQ${params.filterQ}.fastq.gz
    """
}

process cutadapt {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/trimmed_primers_FASTQ", pattern: '*.fastq.gz', mode: params.publish_dir_mode
    publishDir "$params.outdir/cutadapt_summary", pattern: '*.report', mode: params.publish_dir_mode
    cpus params.cutadapt_cpu

    input:
    tuple val(sampleID), path(sampleFASTQ)
    val forward_primer
    val reverse_primer

    output:
    tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq
    path "*.report", emit: cutadapt_summary
    path "cutadapt_summary_${sampleID}.tsv", emit: summary_tocollect
    path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq_files

    script:
    """
    cutadapt -g "${forward_primer}...${reverse_primer}" \
        ${sampleFASTQ} \
        -o ${sampleID}.trimmed.fastq.gz \
        -j ${task.cpus} --trimmed-only --revcomp -e 0.1 \
        --json ${sampleID}.cutadapt.report

    input_read=`jq -r '.read_counts | .input' ${sampleID}.cutadapt.report`
    demux_read=`jq -r '.read_counts | .output' ${sampleID}.cutadapt.report`
    echo -e "sample\tinput_reads\tdemuxed_reads" > cutadapt_summary_${sampleID}.tsv
    echo -e "${sampleID}\t\$input_read\t\$demux_read" >> cutadapt_summary_${sampleID}.tsv
    """
}

// Length filter for the Emu read path. BLCA/NB get their length window enforced
// inside `qiime dada2 denoise-ccs` (--p-min-len/--p-max-len); Emu skips DADA2, so
// without this it would have no length bound. Applying the SAME min_len/max_len
// here keeps all three classifiers on identically filtered reads.
process length_filter {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"
    label 'cpu8'
    publishDir "$params.outdir/emu_length_filtered_FASTQ", pattern: '*.lenfilt.fastq.gz', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(sampleFASTQ)
    val min_len
    val max_len

    output:
    tuple val(sampleID), path("${sampleID}.lenfilt.fastq.gz"), emit: filtered

    script:
    """
    seqkit seq -j ${task.cpus} -m ${min_len} -M ${max_len} ${sampleFASTQ} \
        -o ${sampleID}.lenfilt.fastq.gz
    """
}

process QC_fastq_post_trim {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"
    label 'cpu8'
    publishDir "$params.outdir/filtered_input_FASTQ", pattern: '*post_trim.tsv', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    path "${sampleID}.seqkit.readstats.post_trim.tsv", emit: all_seqkit_stats

    script:
    """
    seqkit fx2tab -j $task.cpus -q --gc -l -H -n -i $sampleFASTQ |\
        csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.readstats.post_trim.tsv
    """
}

process collect_QC {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"
    publishDir "$params.outdir/results/reads_QC", mode: params.publish_dir_mode
    label 'cpu8'

    input:
    path "*"
    path "*"
    path "*"
    path "*"

    output:
    path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
    path "all_samples_seqkit.readstats.post_trim.tsv", emit: all_samples_readstats_post_trim
    path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
    path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
    path "seqkit.summarised_stats.group_by_samples.pretty.tsv"
    path "all_samples_cutadapt_stats.tsv", emit: cutadapt_summary

    script:
    """
    csvtk concat -t -C '%' *.seqkit.readstats.tsv > all_samples_seqkit.readstats.tsv
    csvtk concat -t -C '%' *.seqkit.readstats.post_trim.tsv > all_samples_seqkit.readstats.post_trim.tsv
    csvtk concat -t -C '%' *.seqkit.summarystats.tsv > all_samples_seqkit.summarystats.tsv
    csvtk concat -t cutadapt_summary*.tsv > all_samples_cutadapt_stats.tsv
    # Summary read_qual for each sample
    csvtk summary -t -C '%' -g sample -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median all_samples_seqkit.readstats.tsv > seqkit.summarised_stats.group_by_samples.tsv
    csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv > seqkit.summarised_stats.group_by_samples.pretty.tsv
    """
}

process collect_QC_skip_cutadapt {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"
    publishDir "$params.outdir/results/reads_QC", mode: params.publish_dir_mode
    label 'cpu8'

    input:
    path "*"
    path "*"

    output:
    path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
    path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
    path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
    path "seqkit.summarised_stats.group_by_samples.pretty.tsv"

    script:
    """
    csvtk concat -t -C '%' *.seqkit.readstats.tsv > all_samples_seqkit.readstats.tsv
    csvtk concat -t -C '%' *.seqkit.summarystats.tsv > all_samples_seqkit.summarystats.tsv
    # Summary read_qual for each sample
    csvtk summary -t -C '%' -g sample -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median all_samples_seqkit.readstats.tsv > seqkit.summarised_stats.group_by_samples.tsv
    csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv > seqkit.summarised_stats.group_by_samples.pretty.tsv
    """
}

process prepare_qiime2_manifest {
    label 'cpu_def'
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"

    input: 
    val(samplesFASTQ)
    path(metadata)

    output:
    path "samplefile*.txt", emit: sample_trimmed_file

    """
    echo -e "sample-id\tabsolute-filepath" > samplefile.txt
    echo -e "$samplesFASTQ" | sed -e 's/\\[//g' -e 's/\\]//g' | tr -d '[:space:]' | tr ',' '\\n' | split -l 2 - sample_ind_
    for i in \$(ls sample_ind_*); do sample=\$(cat \${i} | tr '\\n' '\\t' | cut -f1); fastq=\$(basename \$(cat \${i} | tr '\\n' '\\t' | cut -f2)); echo -e "\${sample}\t\${fastq}"; done >> samplefile.txt
    rm -f sample_ind_*
    # If pool column exists, split sample files
    # First make sure metadata file is not empty
    # Check if metadata file exists and is not empty
    if [ ! -s "$metadata" ]; then
        echo "Error: Metadata file is empty or does not exist"
        exit 1
    fi
    poolyes=\$(csvtk headers -t $metadata | grep pool | wc -l)
    if [[ \$poolyes -eq 1 ]]
    then
        csvtk join -t samplefile.txt <(csvtk cut -t -f sample_name,pool $metadata) -f "sample-id;sample_name" | csvtk split -t -f pool
        counter=1
        for i in \$(ls stdin*.tsv); do csvtk cut -t -f -pool \${i} > samplefile\${counter}.txt; rm -f \${i}; let counter++; done
        rm -f samplefile.txt
    fi
    """
}

process prepare_qiime2_manifest_skip_cutadapt {
    label 'cpu_def'
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"

    input: 
    val(samplesFASTQ)
    path(metadata)

    output:
    path "samplefile*.txt", emit: sample_trimmed_file

    """
    echo -e "sample-id\tabsolute-filepath" > samplefile.txt
    echo -e "$samplesFASTQ" | sed -e 's/\\[//g' -e 's/\\]//g' | tr -d '[:space:]' | tr ',' '\\n' | split -l 2 - sample_ind_
    for i in \$(ls sample_ind_*); do sample=\$(cat \${i} | tr '\\n' '\\t' | cut -f1); fastq=\$(basename \$(cat \${i} | tr '\\n' '\\t' | cut -f2)); echo -e "\${sample}\t\${fastq}"; done >> samplefile.txt
    rm -f sample_ind_*
    # If pool column exists, split sample files
    # First make sure metadata file is not empty
    # Check if metadata file exists and is not empty
    if [ ! -s "$metadata" ]; then
        echo "Error: Metadata file is empty or does not exist"
        exit 1
    fi
    poolyes=\$(csvtk headers -t $metadata | grep pool | wc -l)
    if [[ \$poolyes -eq 1 ]]
    then
        csvtk join -t samplefile.txt <(csvtk cut -t -f sample_name,pool $metadata) -f "sample-id;sample_name" | csvtk split -t -f pool
        counter=1
        for i in \$(ls stdin*.tsv); do csvtk cut -t -f -pool \${i} > samplefile\${counter}.txt; rm -f \${i}; let counter++; done
        rm -f samplefile.txt
    fi
    """
}

process import_qiime2 {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/import_qiime", mode: params.publish_dir_mode
    label 'cpu_def'


    input:
    path sample_manifest
    path '*'

    output:
    path 'samples.qza'

    script:
    """
    # Make sure path is absolute 
    export HOME="\$PWD/.qiime_home"   # \$HOME is read-only in the image; redirect all caches here
    mkdir -p "\$HOME"
    awk -v wdir="\$(pwd)/" -F\$'\t' '{if (NR>1){print \$1,wdir\$2} else {print \$0}}' OFS=\$'\t' $sample_manifest > sample_list.txt
    qiime tools import --type 'SampleData[SequencesWithQuality]' \
        --input-path sample_list.txt \
        --output-path samples.qza \
        --input-format SingleEndFastqManifestPhred33V2
    """
}

process demux_summarize {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/summary_demux", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path samples_qza

    output:
    path "samples.demux.summary.qzv"
    path "per-sample-fastq-counts.tsv"

    script:
    """
    export HOME="\$PWD/.qiime_home"   # \$HOME is read-only in the image; redirect all caches here
    mkdir -p "\$HOME"
    qiime demux summarize --i-data $samples_qza \
        --o-visualization samples.demux.summary.qzv

    qiime tools export --input-path samples.demux.summary.qzv \
        --output-path ./
    """
}

process merge_manifest {
    label 'cpu_def'
    publishDir "$params.outdir/results/", mode: params.publish_dir_mode
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"

    input:
    path "*"

    output:
    path "samplefile_merged.txt", emit: merged_manifest

    script:
    """
    # Merge all sample manifests
    echo -e "sample-id\tabsolute-filepath" > samplefile_merged.txt
    for i in \$(ls samplefile*.txt);
    do
        if [ \$i != "samplefile_merged.txt" ];
        then
            awk 'NR>1' \$i >> samplefile_merged.txt
        fi
    done
    """
}

process cutadapt_stats {
    conda (params.enable_conda ? "$projectDir/envs/qiime2-amplicon.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    cpus 1

    input:
    path "cutadapt_summary_*.tsv"

    output:
    path "primer_removal_stats.tsv", emit: stats

    script:
    """
    # merge per-sample cutadapt summaries (sample, input_reads, demuxed_reads)
    awk -F'\\t' 'FNR==1 && seen{next} FNR==1{seen=1} {print}' cutadapt_summary_*.tsv > _all.tsv
    # compute primers removed and % reads with primer (+ a TOTAL row)
    awk -F'\\t' 'BEGIN{OFS="\\t"}
      NR==1{print "sample","input_reads","with_primer","reads_dropped_no_primer","pct_with_primer"; next}
      {inp=\$2+0; kept=\$3+0;
       printf "%s\\t%d\\t%d\\t%d\\t%.2f\\n", \$1, inp, kept, inp-kept, (inp>0?100*kept/inp:0);
       ti+=inp; tk+=kept}
      END{printf "TOTAL\\t%d\\t%d\\t%d\\t%.2f\\n", ti, tk, ti-tk, (ti>0?100*tk/ti:0)}' _all.tsv > primer_removal_stats.tsv
    """
}
