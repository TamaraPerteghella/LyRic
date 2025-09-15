#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import everything from helpers.nf
include { fastqStats } from './processes/fastqStats.nf'

workflow.onComplete {
    println "LyRic workflow finished smoothly."
}

workflow.onError { t ->        // t is a Throwable
    if (t != null) {
        println "LyRic workflow finished with errors: ${t.message}"
    } else {
        println "LyRic workflow finished with errors."
    }
}

// Define the main workflow
workflow {

    def fastq_ch = Channel.fromPath("${params.datadir}/*.fastq.gz")
        .map { file ->
            def base = file.name.replaceAll(/\.fastq\.gz$/, '')
            tuple(base, file)
        }

    // Call the process or workflow from helpers.nf
    fastqStats(fastq_ch)
}