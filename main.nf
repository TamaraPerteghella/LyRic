#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
nextflow.preview.output = true

// Import everything from helpers.nf
include { fastqStats } from './processes/fastqStats.nf'


// Define the main workflow
workflow {
    main:

    def annotation = file("${launchDir}/config/sample_annotations.tsv")

    if (annotation.exists()) {
        println("Using files listed in: ${annotation}")

        fastq_ch = Channel.fromPath("${launchDir}/config/sample_annotations.tsv")
            .splitText()
            .filter { !it.startsWith("absolute_path") }
            .map { line ->
                def (filePath, filename) = line.split('\t')[0..1]
                tuple(filename, file(filePath))
            }
    }
    else {
        println("Expecting files in: ${params.datadir}/")
        fastq_ch = Channel.fromPath("${params.datadir}/*.fastq.gz")
            .map { file ->
                def base = file.name.replaceAll(/\.fastq\.gz$/, '')
                tuple(base, file)
            }
    }


    workflow.onComplete {
        println("LyRic workflow finished smoothly.")
    }

    workflow.onError { t ->
        if (t != null) {
            println("LyRic workflow finished with errors: ${t.message}")
        }
        else {
            println("LyRic workflow finished with errors.")
        }
    }

    // Call the process or workflow from helpers.nf
    fastqStats(fastq_ch)

    statfiles = fastqStats.out.timestamp_fastq.concat(fastqStats.out.aggregate_summary)
    qcfiles = fastqStats.out.qc_ch
    plotfiles = fastqStats.out.read_length_plot
    tmpstatsfile = fastqStats.out.readlength_ch.map { t -> t[1] }.concat(fastqStats.out.readlength_summary_ch)

    publish:
    qc = qcfiles
    stats = statfiles
    stats_tmp = tmpstatsfile
    plots = plotfiles
}

output {
    qc {
        path { fastqc -> "${params.qcdir}/" }
    }
    plots {
        path { plt -> "${params.plotdir}/" }
    }
    stats {
        path { stat -> "${params.statsdir}/" }
    }
    stats_tmp {
        path { stat -> "${params.statsdir}/tmp/" }
    }
}
