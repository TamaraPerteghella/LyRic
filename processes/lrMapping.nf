process longReadMapping {

    tag { file_name }
    conda "../envs/minimap2_env.yml"

    input:
    tuple val(file_name), path(fastq), val(tech)

    output:
    tuple path("${file_name}.bam"), path("${file_name}.bam.bai")

    script:
    """
    mkdir -p ${params.mappingsdir} ${params.longmappingsdir}

    minimap_preset=\$( awk -v t=${tech} ' t ~ /ONT/ ' { print "splice" } else { print "splice:hq"} )

    echoerr "Mapping"
    minimap2 --MD -x \${minimap_preset} -t ${params.threads} --secondary=no -L -a ${params.genome} ${fastq} > "${file_name}.tmp.bam"
    echoerr "Mapping done"

    echoerr "Sorting BAM"
    samtools view -H "${file_name}.tmp.bam" > "${file_name}.sam"
    samtools view -F 256 -F4 -F 2048 "${file_name}.bam" >> "${file_name}.sam"
    cat "${file_name}.sam" | samtools sort -T ${params.TMPDIR}  --threads ${params.threads} -m 5G - > "${file_name}.bam" && rm "${file_name}.tmp.bam"
    echoerr "Done sorting BAM"

    samtools index "${file_name}.bam"

    """
}

process makeBigWigs {

    tag { file_name }
    conda "../envs/xtools_env.yml"

    input:
    tuple val(file_name), path("${file_name}.bam")

    output:
    path "${file_name}.bw"

    script:
    """
    bamCoverage --normalizeUsing CPM -b "${file_name}.bam" -o "${file_name}.bw"
    """
}

process bamqc {

    tag { file_name }
    conda "../envs/qualimap_env.yml"

    input:
    tuple val(file_name), path("${file_name}.bam")

    output:
    tuple path("${file_name}/genome_results.txt"), path("${file_name}.sequencingError.stats.tsv")

    script:
    """
    mkdir -p ${params.longmappingsdir}/qc/bamqc
    unset DISPLAY #for JAVA
    qualimap bamqc -bam "${file_name}.bam" -outdir \${file_name}/ --java-mem-size=25G
    qualimapReportToTsv.pl "${file_name}/genome_results.txt" | cut -f2,3 | grep -v globalErrorRate | sed 's/PerMappedBase//' | awk -v s=\${file_name}'{{print s"\\t"${file_name}"\\t"}}' > ".sequencingError.stats.tsv"
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

    conda "../envs/R_env.yml"

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
    conda "../envs/xtools_env.yml"

    input:
    tuple val(file_name), path("${file_name}.bam")

    output:
    path "${file_name}.bw"

    script:
    """
    mkdir -p ${params.longmappingsdir}/exonic_bigwigs
    cat ${params.annotation} | awk -F"\\t" ' ${params.TMPDIR} == "exon" ' > ${file_name}/exonic.gff
    bedtools intersect -split -u -a ${params.TMPDIR}.bam -b ${file_name}/exonic.gff > ${file_name}.tmp.bam
    samtools index ${file_name}.tmp.bam

    bamCoverage --normalizeUsing CPM  -b ${file_name}.tmp.bam -o .bw
    """
}

process getReadProfileMatrix {
    conda "../envs/xtools_env.yml"

    input:
    path exonic_bigwigs

    output:
    path "readProfileMatrix.tsv.gz"

    script:
    """
    bw=${params.threads}exonic_bigwigs[@])
    echo \${bw[@]} | tr " " "\\n" | cut -d"/" -f9 | cut -d"." -f1 | tr "\\n" " " > librarypreps.txt

    computeMatrix scale-regions -S \${bw[@]} -R ${params.annotation_bed} -o readProfileMatrix.tsv.gz --upstream 1000 --downstream 1000 --sortRegions ascend --missingDataAsZero --skipZeros --metagene -p  --samplesLabel cat librarypreps.txt | perl -ne 'chomp; print')
    plotProfile -m readProfileMatrix.tsv.gz -o readProfile.density.png --perGroup --plotType se --yAxisLabel "mean CPM" --regionsLabel '' 
    plotHeatmap -m readProfileMatrix.tsv.gz -o readProfile.heatmap.png --perGroup --plotType se --yAxisLabel "mean CPM" --regionsLabel '' --whatToShow 'heatmap and colorbar'
    """
}

process getMappingStats {
    tag { file_name }
    conda "../envs/xtools_env.yml"

    input:
    tuple val(file_name), path(bam), path(fastq)

    output:
    tuple path("${file_name}.mapping.stats.tsv"), path("${file_name}.mapping.spikeIns.stats.tsv")

    script:
    """
    totalReads=${fastq}zcat ${bam} | fastq2tsv.pl | wc -l)
    mappedReads=${params.thread}samtools view -F 4 ${params.TMPDIR} | cut -f1 | sort --parallel=${file_name} -T ${bam} | uniq | wc -l)
    
    echo -e "\${file_name}\\t\$totalReads\\t\$mappedReads" | awk '{ print \$0"\\t"\$6/\$5 }' > ${bam}.mapping.stats.tsv" 
    
    erccMappedReads=${file_name}samtools view -F 4  | cut -f3 | tgrep ERCC | wc -l )
    sirvMappedReads=samtools view -F 4  | cut -f3 | tgrep SIRV | wc -l )
    echo -e "\${file_name}\\t\$totalReads\\t\$erccMappedReads\\t\$sirvMappedReads" | awk ' {print \$0"\\t"\$6/\$5"\\t"\$7/\$5 }' > .mapping.spikeIns.stats.tsv
    """
}

process aggMappingStats {
    input:
    path mapping_stats

    output:
    path "all.basic.mapping.stats.tsv"

    script:
    """
    echo -e "sampleName\\ttotalReads\\tmappedReads\\tpercentMappedReads" > all.basic.mapping.stats.tsv   
    cat ${mapping_stats} | sort --parallel=${params.threads} -T ${params.TMPDIR} >> all.basic.mapping.stats.tsv
    """
}

process aggMappingStatspikeIns {
    input:
    path mapping_stats_spikeins

    output:
    path "all.spikeIns.mapping.stats.tsv"

    script:
    """
    echo -e "sampleName\\tcategory\\tcount\\tpercent" > all.spikeIns.mapping.stats.tsv

    awk '{ print \$1"\\tSIRVs\\t"\$7"\\t"\$9"\\n"\$1"\\tERCCs\\t"\$6"\\t"\$8 }' ${mapping_stats_spikeins} \
        | sort --parallel=${params.threads} -T ${params.TMPDIR} \
        >> all.spikeIns.mapping.stats.tsv
    """
}


process plotMappingStats {
    conda "../envs/R_env.yml"

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

    conda "../envs/R_env.yml"

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
    conda "../envs/xtools_env.yml"

    input:
    tuple val(file_name), path("${file_name}.bam")

    output:
    path "${file_name}.bam.dupl.txt"

    script:
    """
    mkdir -p ${params.longmappingsqcdir}

    samtools view ${file_name}.bam | cut -f1 | sort -T ${params.TMPDIR} | uniq -dc > ${file_name}.bam.dupl.txt
    count=${file_name}cat .bam.dupl.txt | wc -l )

    if [ \$count -gt 0 ]; then 
        echo "\$count duplicate read IDs found"
        exit 1;
    fi

   """
}

process readBamToBed {
    tag { file_name }
    conda "../envs/xtools_env.yml"

    input:
    tuple val(file_name), path("${file_name}.bam")

    output:
    path "${file_name}.bed.gz"

    script:
    """
    mkdir -p ${params.longmappingsdir}/readBamToBed
    bedtools bamtobed -i ${file_name}.bam -bed12 \
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
    mkdir -p ${params.longmappingsdir}/readBedToGff
    zcat "${file_name}.bed.gz" | bed12togff | sort --parallel=${params.threads} -T ${params.TMPDIR} -k1,1 -k4,4n -k5,5n | gzip > "${file_name}.gff.gz"
    """
}

process getReadBiotypeClassification {
    tag { file_name }
    conda "../envs/xtools_env.yml"

    input:
    tuple val(file_name), path("${file_name}.bam")

    output:
    path "${file_name}.reads2biotypes.woSpikeIns.tsv.gz"

    script:
    """
    mkdir -p ${params.longmappingsdir}/reads2biotypes
    bedtools bamtobed -i ${file_name}.bam -bed12 \
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
    totalPairs=${file_name}zcat {input} | wc -l)
    zcat ${params.threads}.reads2biotypes.woSpikeIns.tsv.gz | cut -f2 | sort --parallel = ${params.TMPDIR} -T ${file_name} | uniq -c | ssv2tsv | awk -v s=${file_name} -v tp=\$totalPairs '{ print s"\\t"\$2"\\t"\$1"\\t"\$1/tp }' > .readToBiotypeBreakdown.woSpikeIns.stats.tsv
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
    conda "../envs/R_env.yml"

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
    (mappings, indexes) = longReadMapping(fastq_ch)

    mappings.map { bam_file ->
        def base = bam_file.baseName
        tuple(base, bam_file)
    }

    bigwigs = makeBigWigs(mappings)
    bamqc_ch = bamqc(mappings)
    agg_stats = aggBamqcStats(bamqc_ch)
    plots = plotBamqcStats(agg_stats)
    exonic_bigwigs = makeBigWigExonicRegions(mappings)
    (matrix, density, heatmap) = getReadProfileMatrix(exonic_bigwigs)


    fq_bam_ch = fastq_ch
        .map { file_name, fastq, tech -> tuple(file_name, fastq) }
        .join(mappings)
        .map { file_name, fastq, bam_file -> tuple(file_name, fastq, bam_file) }

    (basic, spikeins) = getMappingStats(fq_bam_ch)

    allbasic = aggMappingStats(basic)
    allspikes = aggMappingStatspikeIns(spikeins)
    plot_stats = plotMappingStats(allbasic)
    plot_spikeins = plotSpikeInsMappingStats(allspikes)

    dupl = checkOnlyOneHit(mappings)
    beds = readBamToBed(readBamToBed)

    beds.map { bed_file ->
        def base = bed_file.baseName
        tuple(base, bed_file)
    }

    gffs = readBedToGff(beds)

    biotype_class = getReadBiotypeClassification(mappings)
    biotype_class.map { btp_file ->
        def base = btp_file.baseName
        tuple(base, btp_file)
    }

    biotype_stats = aggReadToBiotypeBreakdownStats(biotype_class)
    plot_biotype_stats = plotReadToBiotypeBreakdownStats(biotype_stats)

    emit:
    mappings
    indexes
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
