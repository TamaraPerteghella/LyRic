// ======================
// Define stats files as a single channel
// ======================

all_stats_ch = Channel.value([
    sampleAnnotationFile     : file("${launchDir}/config/sample_annotations.tsv"),
    allFastqTimeStamps       : file("${params.statsdir}/all.fastq.timestamps.tsv"),
    allReadLengths           : file("${params.statsdir}/all.readlength.summary.tsv"),
    allBasicMappingStats     : file("${params.statsdir}/all.basic.mapping.stats.tsv"),
    allHissStats             : file("${params.statsdir}/all.HiSS.stats.tsv"),
    allMergedStats           : file("${params.statsdir}/all.min${params.minReadSupport}reads.merged.stats.tsv"),
    allMatureRnaLengthStats  : file("${params.statsdir}/all.min${params.minReadSupport}reads.matureRNALengthSummary.stats.tsv"),
    allTmergeVsSirvStats     : file("${params.statsdir}/all.HiSS.tmerge.min${params.minReadSupport}reads.vs.SIRVs.stats.tsv"),
    allCagePolyASupportStats : file("${params.statsdir}/all.min${params.minReadSupport}reads.splicing_status-all.cagePolyASupport.stats.tsv"),
    allNovelLociStats        : file("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-all.novelLoci.stats.tsv"),
    allNovelFlLociStats      : file("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-cagePolyASupported.novelLoci.stats.tsv"),
    allNovelLociQcStats      : file("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-all.novelLoci.qc.stats.tsv"),
    allNovelFlLociQcStats    : file("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-cagePolyASupported.novelLoci.qc.stats.tsv"),
    allNtCoverageStats       : file("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-all.vs.ntCoverageByGenomePartition.stats.tsv")
])

// ======================
// Filter sample annotation for a given project
// ======================
process filterSampleAnnot {
    tag { prj }

    input:
    path sample_annotation
    val prj

    output:
    path "${prj}_samples.tsv"

    script:
    """
    grep -v '^absolute_path' ${sample_annotation} | awk -F'\\t' '\$4=="${prj}" {print \$2}' > ${prj}_samples.tsv
    """
}

// ======================
// Make HTML summary dashboard
// ======================
process makeHtmlSummaryDashboard {
    tag { prj }
    conda "../envs/R_env.yml"

    input:
    val stats_map from all_stats_ch
    val prj

    output:
    path("${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.html")
    path("${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.tsv")
    path("${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.index.tmp.html")

    script:
    """
    makeHtmlDashboard.r ${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.html \\
        ${stats_map.sampleAnnotationFile} ${stats_map.allFastqTimeStamps} ${stats_map.allReadLengths} \\
        ${stats_map.allBasicMappingStats} ${stats_map.allHissStats} ${stats_map.allMergedStats} \\
        ${stats_map.allMatureRnaLengthStats} ${stats_map.allTmergeVsSirvStats} ${stats_map.allCagePolyASupportStats} \\
        ${stats_map.allNovelLociStats} ${stats_map.allNovelFlLociStats} ${stats_map.allNovelLociQcStats} \\
        ${stats_map.allNovelFlLociQcStats} ${stats_map.allNtCoverageStats}

    htmlBn="${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}"
    tsvBn="${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}"

    echo "<li> <b>${prj}</b> ${params.indexEntryPart}: <a href='\$htmlBn'>HTML</a> / <a href='\$tsvBn'>TSV</a></li>" \
        > "${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.index.tmp.html"
    """
}

// ======================
// Make HTML index page
// ======================
process makeHtmlSummaryDashboardIndex {
    input:
    path index_entries
    val prj

    output:
    path "output/html/index.html"

    script:
    """
    #mkdir -p output/html
    printf '<!DOCTYPE html>\\n' > output/html/index.html
    printf "<html>\\n<head><title>Summary statistics tables for ${prj} project</title></head><body>\\n" >> output/html/index.html
    printf "<h1>Summary statistics tables for ${prj} project</h1>\\n<ul>\\n" >> output/html/index.html

    cat ${index_entries} | sort >> output/html/index.html

    printf '</ul>' >> output/html/index.html
    printf "<br>Produced with <a href='https://github.com/julienlag/LyRic/'>LyRic</a>.<br>" >> output/html/index.html
    date=$(date)
    printf "<br>(Last updated \$date)" >> output/html/index.html
    printf '</body></html>' >> output/html/index.html
    """
}

// ======================
// Workflow
// ======================
workflow htmlTables {

    if (params.split_html) {

        // Read sample annotation and extract project info
        projects_ch = Channel.fromPath("${launchDir}/config/sample_annotations.tsv")
            .splitText()
            .filter { !it.startsWith("absolute_path") }
            .map { line ->
                def (filePath, prj) = line.split('\t')[1,3]
                tuple(filePath, prj)
            }

        // Get unique projects
        projects_ch
            .map { filePath, prj -> prj }
            .distinct()
            .subscribe { prj ->
                // Filter samples per project
                filterSampleAnnot(file("${launchDir}/config/sample_annotations.tsv"), prj)
                makeHtmlSummaryDashboard(all_stats_ch, prj)
            }

    } else {
        // Single merged report for all projects
        projects_ch = Channel.fromPath("${launchDir}/config/sample_annotations.tsv")
        projects_ch.toPath("samples.tsv")
        makeHtmlSummaryDashboard(all_stats_ch, "all_projects")
    }
}
