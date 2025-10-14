process longReadMapping {
    tag { file_name }

    input:
    tuple val(file_name), path(fastq), val(genome), val(tech)

    output:
    tuple path("${file_name}.bam"), path("${file_name}.bam.bai")

    script:
    """
    if [[ ${tech} == "ONT" ]]; then 
        minimap_preset="splice" 
    else 
        minimap_preset="splice:hq" 
    fi

    echoerr "Mapping"
    minimap2 --MD -x \${minimap_preset} -t ${task.cpus} --secondary=no -L -a ${params.genomes_folder}/${genome}.fa.gz ${fastq} > "${file_name}.tmp.bam"
    echoerr "Mapping done"

    echoerr "Sorting BAM"
    samtools view -H "${file_name}.tmp.bam" > "${file_name}.sam"
    samtools view -F 256 -F4 -F 2048 "${file_name}.tmp.bam" >> "${file_name}.sam"
    samtools sort -T ${params.TMPDIR} --threads ${params.threads} -m 2G "${file_name}.sam" > "${file_name}.bam" && rm "${file_name}.tmp.bam"
    echoerr "Done sorting BAM"

    samtools index "${file_name}.bam"

    """
}

process makeBigWigs {
    tag { file_name }

    input:
    tuple val(file_name), path(bam_file), path(bai_file)

    output:
    path "${file_name}.bw"

    script:
    """
    bamCoverage --normalizeUsing CPM -b ${bam_file} -o "${file_name}.bw"
    """
}

process bamqc {
    tag { file_name }

    input:
    tuple val(file_name), path(bam_file), path(bai_file)

    output:
    tuple path("${file_name}/genome_results.txt"), path("${file_name}.sequencingError.stats.tsv")

    script:
    """
    unset DISPLAY #for JAVA
    qualimap bamqc -bam ${bam_file} -outdir ${file_name}/ --java-mem-size=25G
    qualimapReportToTsv.pl "${file_name}/genome_results.txt" | cut -f2,3 | grep -v globalErrorRate | sed 's/PerMappedBase//' | awk -v s=${file_name} '{ print s"\\t"\$1"\\t"\$2} ' > "${file_name}.sequencingError.stats.tsv"
    """
}

process aggBamqcStats {
    input:
    path sequencingError_file

    output:
    path "all.sequencingError.stats.tsv"

    script:
    """
    echo -e "sample_name\\terrorCategory\\terrorRate" > all.sequencingError.stats.tsv
    cat ${sequencingError_file} | sort --parallel=${params.threads} >> all.sequencingError.stats.tsv
    """
}

process plotBamqcStats {
    label 'rplots'

    input:
    path seq_error_stats

    output:
    tuple path("sequencingError.allErrors.stats.pdf"), path("sequencingError.deletionsOnly.stats.pdf")

    script:
    """
    plotBamqcStats.R ${seq_error_stats} sequencingError.allErrors.stats.pdf sequencingError.deletionsOnly.stats.pdf
    """
}



process makeBigWigExonicRegions {
    tag { file_name }

    input:
    tuple val(file_name), path(bam_file), path(bai_file)

    output:
    path "${file_name}.bw"

    script:
    """
    awk -F"\\t" '\$3 == "exon" ' ${params.annotation} > ${file_name}.exonic.gff
    bedtools intersect -split -u -a ${bam_file} -b ${file_name}.exonic.gff > ${file_name}.tmp.bam
    samtools index ${file_name}.tmp.bam

    bamCoverage --normalizeUsing CPM -b ${file_name}.tmp.bam -o ${file_name}.bw
    """
}

process getReadProfileMatrix {
    input:
    path exonic_bigwigs

    output:
    tuple path ("readProfileMatrix.tsv.gz"), path("readProfile.density.png"), path("readProfile.heatmap.png")

    script:
    """
    bw=( ${exonic_bigwigs} )
    echo \${bw[@]} | tr " " "\\n" | cut -d"/" -f9 | cut -d"." -f1 | tr "\\n" " " > librarypreps.txt

    computeMatrix scale-regions -S \${bw[@]} -R ${params.annotation} -o readProfileMatrix.tsv.gz --upstream 1000 --downstream 1000 --sortRegions ascend --missingDataAsZero --skipZeros --metagene -p ${params.threads} --samplesLabel \$(cat librarypreps.txt | perl -ne 'chomp; print')
    plotProfile -m readProfileMatrix.tsv.gz -o readProfile.density.png --perGroup --plotType se --yAxisLabel "mean CPM" --regionsLabel '' 
    plotHeatmap -m readProfileMatrix.tsv.gz -o readProfile.heatmap.png --perGroup --plotType se --yAxisLabel "mean CPM" --regionsLabel '' --whatToShow 'heatmap and colorbar'
    """
}

process getMappingStats {
    tag { file_name }

    input:
    tuple val(file_name), path(fastq), path(bam)

    output:
    tuple path("${file_name}.mapping.stats.tsv"), path("${file_name}.mapping.spikeIns.stats.tsv")

    script:
    """
    totalReads=\$(zcat ${fastq} | fastq2tsv.pl | wc -l)
    mappedReads=\$(samtools view -F 4 ${bam} | cut -f1 | sort --parallel=${params.threads} -T ${params.TMPDIR} | uniq | wc -l)
    
    echo -e "${file_name}\\t\$totalReads\\t\$mappedReads" | awk '{ print \$0"\\t"\$3/\$2 }' > "${file_name}.mapping.stats.tsv"
    
    erccMappedReads=\$(samtools view -F 4 ${bam} | cut -f3 | tgrep ERCC | wc -l )
    sirvMappedReads=\$(samtools view -F 4 ${bam} | cut -f3 | tgrep SIRV | wc -l )
    echo -e "${file_name}\\t\$totalReads\\t\$erccMappedReads\\t\$sirvMappedReads" | awk ' {print \$0"\\t"\$3/\$2"\\t"\$4/\$2 }' > "${file_name}.mapping.spikeIns.stats.tsv"
    """
}

process aggMappingStats {
    input:
    path mapping_stats

    output:
    path "all.basic.mapping.stats.tsv"

    script:
    """
    echo -e "sample_name\\ttotalReads\\tmappedReads\\tpercentMappedReads" > all.basic.mapping.stats.tsv   
    cat ${mapping_stats} | sort --parallel=${params.threads} -T ${params.TMPDIR} >> all.basic.mapping.stats.tsv
    """
}

process aggMappingStatspikeIns {
    tag "Aggregate mapping Stat for Spike-Ins"

    input:
    path mapping_stats_spikeins

    output:
    path "all.spikeIns.mapping.stats.tsv"

    script:
    """
    echo -e "sample_name\\tcategory\\tcount\\tpercent" > all.spikeIns.mapping.stats.tsv

    awk '{ print \$1"\\tSIRVs\\t"\$7"\\t"\$9"\\n"\$1"\\tERCCs\\t"\$6"\\t"\$8 }' ${mapping_stats_spikeins} \
        | sort --parallel=${params.threads} -T ${params.TMPDIR} >> all.spikeIns.mapping.stats.tsv
    """
}


process plotMappingStats {
    label 'rplots'

    input:
    path basic_stats

    output:
    path "lrMapping.basic.stats"

    script:
    """
    plotMappingStats.R ${basic_stats} lrMapping.basic.stats
    """
}

process plotSpikeInsMappingStats {
    label 'rplots'

    input:
    path spikeins_stats

    output:
    path "lrMapping.spikeIns.stats"

    script:
    """
    plotSpikeInsMappingStats.R ${spikeins_stats} lrMapping.spikeIns.stats
    """
}

process checkOnlyOneHit {
    tag { file_name }

    input:
    tuple val(file_name), path(bam_file), path(bai_file)

    output:
    path "${file_name}.bam.dupl.txt"

    script:
    """
    samtools view ${bam_file} | cut -f1 | sort --parallel=${params.threads} -T ${params.TMPDIR} | uniq -dc > ${file_name}.bam.dupl.txt
    count=\$(cat ${file_name}.bam.dupl.txt | wc -l )

    if [ \$count -gt 0 ]; then 
        echo "\$count duplicate read IDs found"
        exit 1;
    fi

   """
}

process readBamToBed {
    tag { file_name }

    input:
    tuple val(file_name), path(bam_file), path(bai_file)

    output:
    path "${file_name}.bed.gz"

    script:
    """
    bedtools bamtobed -i ${bam_file} -bed12 \
      | perl -ne '\$line=\$_; @line=split("\\t", \$line); @blockSizes=split(",", \$line[10]); \$allExonsOK=1; foreach \$block (@blockSizes){ if (\$block<2){ \$allExonsOK=0; last; } }; if (\$allExonsOK==1){ print \$line }' \
      | sort --parallel=${params.threads} -T ${params.TMPDIR} -k1,1 -k2,2n -k3,3n \
      | gzip > ${file_name}.bed.gz
    """
}

process readBedToGff {
    input:
    tuple val(file_name), path("${file_name}.bed.gz")

    output:
    path "${file_name}.gff.gz"

    script:
    """
    zcat "${file_name}.bed.gz" | bed12togff | sort --parallel=${params.threads} -T ${params.TMPDIR} -k1,1 -k4,4n -k5,5n | gzip > "${file_name}.gff.gz"
    """
}

process getReadBiotypeClassification {
    tag { file_name }

    input:
    tuple val(file_name), path(bam_file), path(bai_file)

    output:
    path "${file_name}.reads2biotypes.woSpikeIns.tsv.gz"

    script:
    """
    bedtools bamtobed -i ${bam_file} -bed12 \
      | bedtools intersect -split -wao -bed -a - -b ${params.annotation} | grep -v '^ERCC' | grep -v '^SIRV' \
      | perl -lane '\$gid="NA"; \$gt="nonExonic"; if(/gene_id "(\\S+)";/){\$gid=\$1} if(/gene_type "(\\S+)";/){\$gt=\$1} print "\$F[3]\\t\$gid\\t\$gt\\t\$F[-1]"' \
      | cut -f1,3 | sort --parallel=${params.threads} -T ${params.TMPDIR} | uniq | gzip > ${file_name}.reads2biotypes.woSpikeIns.tsv.gz
    """
}

process getReadToBiotypeBreakdownStats {
    input:
    tuple val(file_name), path("${file_name}.reads2biotypes.woSpikeIns.tsv.gz")

    output:
    path "${file_name}.readToBiotypeBreakdown.woSpikeIns.stats.tsv"

    script:
    """
    totalPairs=\$( zcat ${file_name}.reads2biotypes.woSpikeIns.tsv.gz | wc -l)
    zcat ${file_name}.reads2biotypes.woSpikeIns.tsv.gz | cut -f2 | sort --parallel=${params.threads} -T ${params.TMPDIR} | uniq -c | ssv2tsv | awk -v s=${file_name} -v tp=\$totalPairs '{ print s"\\t"\$2"\\t"\$1"\\t"\$1/tp }' > ${file_name}.readToBiotypeBreakdown.woSpikeIns.stats.tsv
    """
}

process aggReadToBiotypeBreakdownStats {
    input:
    path stats_biotypes

    output:
    path "all.readToBiotypeBreakdown.woSpikeIns.stats.tsv"

    script:
    """
    echo -e "sample_name\\tbiotype\\treadOverlapsCount\\treadOverlapsPercent" > all.readToBiotypeBreakdown.woSpikeIns.stats.tsv
    cat ${stats_biotypes} | sort --parallel=${params.threads} -T ${params.TMPDIR} >> all.readToBiotypeBreakdown.woSpikeIns.stats.tsv
    """
}

process plotReadToBiotypeBreakdownStats {
    label 'rplots'
    
    input:
    path stats_biotypes_out

    output:
    path "readToBiotypeBreakdown.stats.woSpikeIns.pdf"

    script:
    """
    plotReadToBiotypeBreakdownStats.R ${stats_biotypes_out} readToBiotypeBreakdown.stats.woSpikeIns.pdf
    """
}

workflow lrMapping {
    take:
    fastq_ch

    main:
    mappings_output = longReadMapping(fastq_ch)

    mappings = mappings_output.map { bam, bai ->
        def file_name = bam.baseName
        tuple(file_name, bam, bai)
    }

    bigwigs = makeBigWigs(mappings)
    bamqc_ch = bamqc(mappings)
    agg_stats = aggBamqcStats(bamqc_ch.map{ _gen, seq -> seq }.collect())
    plots = plotBamqcStats(agg_stats)
    exonic_bigwigs = makeBigWigExonicRegions(mappings)
    (matrix, density, heatmap) = getReadProfileMatrix(exonic_bigwigs.collect())


    fq_bam_ch = fastq_ch
        .map { file_name, fastq, _genome, _tech -> tuple(file_name, fastq) }
        .join(mappings)
        .map { file_name, fastq, bam_file, _bai_file -> tuple(file_name, fastq, bam_file) }

    mappings_stats = getMappingStats(fq_bam_ch)

    allbasic = aggMappingStats(mappings_stats.map{ basic, _spikeins -> basic }.collect())
    allspikes = aggMappingStatspikeIns(mappings_stats.map{ _basic, spikeins -> spikeins }.collect())
    plot_stats = plotMappingStats(allbasic)
    plot_spikeins = plotSpikeInsMappingStats(allspikes)

    dupl = checkOnlyOneHit(mappings)
    beds = readBamToBed(mappings)

    gffs = readBedToGff(beds.map { bed_file -> def base = bed_file.baseName; tuple(base, bed_file) })

    biotype_class = getReadBiotypeClassification(mappings)

    biotype_stats = aggReadToBiotypeBreakdownStats(biotype_class.collect())
    plot_biotype_stats = plotReadToBiotypeBreakdownStats(biotype_stats)

    emit:
    mappings
    bigwigs
    bamqc_ch
    agg_stats
    plots
    exonic_bigwigs
    matrix
    density
    heatmap
    allbasic
    allspikes
    plot_stats
    plot_spikeins
    dupl
    beds
    gffs
    biotype_class
    biotype_stats
    plot_biotype_stats
}
