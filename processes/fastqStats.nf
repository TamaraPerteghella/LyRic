// Check for duplicate IDs in fastq
process basicFASTQqc {

    tag { file_name }

    input:
    tuple val(file_name), path(fastq), val(tech)

    output:
    path "${file_name}.dupl.txt"

    script:
    """
    mkdir -p ${params.qcdir}

    # check that there are no read ID duplicates
    zcat ${fastq} | fastq2tsv.pl | awk '{print \$1}' | sort --parallel=${params.threads} -T ${params.TMPDIR} | uniq -dc > ${file_name}.dupl.txt

    count=\$(cat ${file_name}.dupl.txt | wc -l)
    if [ \$count -gt 0 ]; then
        echo "\$count duplicate read IDs found"
        mv ${file_name}.dupl.txt ${file_name}.dupl.txt.tmp
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
    path "all.fastq.timestamps.tsv"

    script:
    """
    mkdir -p ${params.statsdir}
    echo -e "sample_name\tFASTQ_modified" > all.fastq.timestamps.tsv

    for t in ${dupl_files}; do
        name=\$(basename \$t .dupl.txt)
        moddate=\$(date -r \$t +%F)
        echo -e "\$name\\t\$moddate" >> all.fastq.timestamps.tsv
    done
    """
}


// get read lengths for all FASTQ files:
process getReadLengthSummary {
    tag { file_name }

    input:
    tuple val(file_name), path(fastq), val(tech)

    output:
    tuple val(file_name), path("${file_name}.readlength.tsv.gz")
    path "${file_name}.readlengthSummary.tsv"

    script:
    """
    echo -e "sample_name\tlength" > ${file_name}.readlength.tsv
    zcat ${fastq} | fastq2tsv.pl | perl -F"\\t" -slane '\$F[0]=~s/^(\\S+).*/\$1/; print join("\\t", @F)' | awk -v s=${file_name} '{print s, length(\$2)}' OFS="\\t" >> ${file_name}.readlength.tsv
    gzip ${file_name}.readlength.tsv

    getReadLengthSummary.R ${file_name}.readlength.tsv.gz ${file_name}.readlengthSummary.tsv
    """
}


process aggReadLengthSummary {
    tag "Aggregate read length summaries"

    input:
    path readlength_files

    output:
    path "all.readlength.summary.tsv"

    script:
    """
    head -n1 ${readlength_files[0]} > all.readlength.summary.tsv
    tail -q -n+2 ${readlength_files} | sort --parallel=${params.threads} >> all.readlength.summary.tsv
    """
}

// plot histograms with R:
process plotReadLength {

    tag { file_name }


    input:
    tuple val(file_name), path(readlength_file)

    output:
    path "${file_name}_readLength.stats.pdf"

    script:
    """
    mkdir -p ${params.plotdir}
    plotReadLength.R ${readlength_file} ${file_name}_readLength.stats.pdf
    """
}


workflow fastqStats {
    take:
    fastq_ch

    main:
    qc_ch = basicFASTQqc(fastq_ch)
    timestamp_fastq = fastqTimestamps(qc_ch)

    (readlength_ch, readlength_summary_ch) = getReadLengthSummary(fastq_ch)
    aggregate_summary = aggReadLengthSummary(readlength_summary_ch)

    read_length_plot = plotReadLength(readlength_ch)

    emit:
    qc_ch
    readlength_ch
    readlength_summary_ch
    aggregate_summary
    read_length_plot
    timestamp_fastq
}
