// Check for duplicate IDs in fastq
process basicFASTQqc {

    tag { file_name }

    input:
    tuple val(file_name), path(fastq)

    output:
    path "${params.qcdir}/${file_name}.dupl.txt"

    script:
    """
    mkdir -p ${params.qcdir}

    # check that there are no read ID duplicates
    zcat ${fastq} | fastq2tsv.pl | awk '{print \$1}' | sort --parallel=${params.threads} -T \$TMPDIR | uniq -dc > ${params.qcdir}/${file_name}.dupl.txt

    count=\$(cat ${params.qcdir}/${file_name}.dupl.txt | wc -l)
    if [ \$count -gt 0 ]; then
        echo "\$count duplicate read IDs found"
        mv ${params.qcdir}/${file_name}.dupl.txt ${params.qcdir}/${file_name}.dupl.txt.tmp
        exit 1
    fi
    """
}


// Generate timestamps TSV
process fastqTimestamps {

    tag "Collect FASTQ modification timestamps"

    input:
    path dupl_files

    output:
    path "${params.statsdir}/all.fastq.timestamps.tsv"

    script:
    """
    mkdir -p ${params.statsdir}
    echo -e "sample_name\tFASTQ_modified" > ${params.statsdir}/all.fastq.timestamps.tsv

    for t in ${dupl_files}; do
        name=\$(basename \$t .dupl.txt)
        moddate=\$(date -r \$t +%F)
        echo -e "\$name\\t\$moddate" >> ${params.statsdir}/all.fastq.timestamps.tsv
    done
    """
}


// get read lengths for all FASTQ files:
process getReadLengthSummary {
    tag { file_name }
    conda "../envs/R_env.yml"

    input:
    tuple val(file_name), path(fastq)

    output:
    tuple val(file_name), path("${params.statsdir}/tmp/${file_name}.readlength.tsv.gz")
    path "${params.statsdir}/tmp/${file_name}.readlengthSummary.tsv"

    script:
    """
    echo -e "sample_name\tlength" > ${params.statsdir}/tmp/${file_name}.readlength.tsv
    zcat ${fastq} | fastq2tsv.pl | perl -F"\\t" -slane '\$F[0]=~s/^(\\S+).*/\$1/; print join("\\t", @F)' | awk -v s=${file_name} '{print s, length(\$2)}' OFS="\\t" >> ${params.statsdir}/tmp/${file_name}.readlength.tsv
    gzip ${params.statsdir}/tmp/${file_name}.readlength.tsv

    Rscript getReadLengthSummary.R ${params.statsdir}/tmp/${file_name}.readlength.tsv.gz ${params.statsdir}/tmp/${file_name}.readlengthSummary.tsv
    """
}


process aggReadLengthSummary {
    tag "Aggregate read length summaries"

    input:
    path readlength_files

    output:
    path "${params.statsdir}/all.readlength.summary.tsv"

    script:
    """
    head -n1 ${readlength_files[0]} > ${params.statsdir}/all.readlength.summary.tsv
    tail -q -n+2 ${readlength_files} | sort --parallel=${params.threads} >> ${params.statsdir}/all.readlength.summary.tsv
    """
}

// plot histograms with R:
process plotReadLength {

    tag { file_name }
    conda "../envs/R_env.yml"

    input:
    tuple val(file_name), path(readlength_file)

    output:
    path "${params.plotdir}/readLength.stats/${file_name}_readLength.stats.pdf"

    script:
    """
    mkdir -p ${params.plotdir}/readLength.stats/
    Rscript plotReadLength.R ${readlength_file} ${params.plotdir}/readLength.stats/${file_name}_readLength.stats.pdf
    """
}


workflow fastqStats {
    take:
    fastq_ch

    main:
    def qc_ch = basicFASTQqc(fastq_ch)
    fastqTimestamps(qc_ch)

    def (readlength_ch, readlength_summary_ch) = getReadLengthSummary(fastq_ch)

    aggReadLengthSummary(readlength_summary_ch)
    plotReadLength(readlength_ch)
}
