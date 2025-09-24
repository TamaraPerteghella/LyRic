// Define channels for each stats file
sampleAnnotationFile_ch     = Channel.fromPath("${launchDir}/config/sample_annotations.tsv")
allFastqTimeStamps_ch       = Channel.fromPath("${params.statsdir}/all.fastq.timestamps.tsv")
allReadLengths_ch           = Channel.fromPath("${params.statsdir}/all.readlength.summary.tsv")
allBasicMappingStats_ch     = Channel.fromPath("${params.statsdir}/all.basic.mapping.stats.tsv")
allHissStats_ch             = Channel.fromPath("${params.statsdir}/all.HiSS.stats.tsv")
allMergedStats_ch           = Channel.fromPath("${params.statsdir}/all.min${params.minReadSupport}reads.merged.stats.tsv")
allMatureRnaLengthStats_ch  = Channel.fromPath("${params.statsdir}/all.min${params.minReadSupport}reads.matureRNALengthSummary.stats.tsv")
allTmergeVsSirvStats_ch     = Channel.fromPath("${params.statsdir}/all.HiSS.tmerge.min${params.minReadSupport}reads.vs.SIRVs.stats.tsv")
allCagePolyASupportStats_ch = Channel.fromPath("${params.statsdir}/all.min${params.minReadSupport}reads.splicing_status-all.cagePolyASupport.stats.tsv")
allNovelLociStats_ch        = Channel.fromPath("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-all.novelLoci.stats.tsv")
allNovelFlLociStats_ch      = Channel.fromPath("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-cagePolyASupported.novelLoci.stats.tsv")
allNovelLociQcStats_ch      = Channel.fromPath("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-all.novelLoci.qc.stats.tsv")
allNovelFlLociQcStats_ch    = Channel.fromPath("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-cagePolyASupported.novelLoci.qc.stats.tsv")
allNtCoverageStats_ch       = Channel.fromPath("${params.statsdir}/all.tmerge.min${params.minReadSupport}reads.endSupport-all.vs.ntCoverageByGenomePartition.stats.tsv")

all_stats_ch = Channel.of(
    tuple(
        sampleAnnotationFile_ch,
        allFastqTimeStamps_ch,
        allReadLengths_ch,
        allBasicMappingStats_ch,
        allHissStats_ch,
        allMergedStats_ch,
        allMatureRnaLengthStats_ch,
        allTmergeVsSirvStats_ch,
        allCagePolyASupportStats_ch,
        allNovelLociStats_ch,
        allNovelFlLociStats_ch,
        allNovelLociQcStats_ch,
        allNovelFlLociQcStats_ch,
        allNtCoverageStats_ch
    )
)

process filterSampleAnnot{
    
    input:
    val(filename), val(prj)
    
    output:
    path "${prj}_samples.csv" 

    script:
    """  
    for f in ${filename};
    do
        echo \${f} >> \${prj}_samples.tsv
    done
    """        
}

process makeHtmlSummaryDashboard{
    conda "../envs/R_env.yml"

    input:
    tuple path(all_stats_ch), val(prj)

    output:
    tuple path("${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.html"), path("${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.tsv"), path("${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.index.tmp.html")

    script:
    """
    makeHtmlDashboard.r ${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.html ${all_stats_ch}
    htmlBn="${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}"
    tsvBn="${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}"

    echo "<li> <b>${prj}</b> ${params.indexEntryPart}: <a href='\$htmlBn'>HTML</a> / <a href='\$tsvBn'>TSV</a></li>" > "${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.index.tmp.html
    """
}

process makeHtmlSummaryDashboardIndex{
    input:
    tuple path("${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.index.tmp.html"), val(prj)

    output:
    "output/html/index.html"

    script:
    """
    # print start of html:
    printf '<!DOCTYPE html\\>' > output/html/index.html
    printf "<html>
    <head>
    <title>Summary statistics tables for {config[PROJECT_NAME]} project</title>
    </head><body>
    <h1>Summary statistics tables for {config[PROJECT_NAME]} project</h1>
    <ul> " >> output/html/index.html

    # print list:
    cat ${params.outputdir}/html/summary_table_min${params.minReadSupport}reads_${prj}.index.tmp.html | sort >> output/html/index.html

    # print end of html:
    printf '</ul>'  >> output/html/index.html

    printf "<br>Produced with <a href="https://github.com/julienlag/LyRic/">LyRic</a>.<br>" >> output/html/index.html
    date=\$(date)
    printf "<br>(Last updated \$date)"  >> output/html/index.html
    printf '</body>
    </html> ' >> output/html/index.html
    """
}

workflow htmlTables {

    if ( ${params.split_html} )
    {
        projects = Channel.fromPath("${launchDir}/config/sample_annotations.tsv")
            .splitText()
            .filter { !it.startsWith("absolute_path") }
            .map { line ->
                def (filePath, projects) = line.split('\t')[1,3]
                tuple(filePath, projects)
            }
        
        for (prj in projects.projects.uniq)
        {
            def tmp_prj = projects.filter{ projects.projects == prj}
            filterSampleAnnot(tmp_prj.filename, prj)
            makeHtmlSummaryDashboard(all_stats_ch, prj)
        }
    }
    else 
    {
        projects.toPath("samples.tsv")
        makeHtmlSummaryDashboard(all_stats_ch, "all_projects")
    }
}