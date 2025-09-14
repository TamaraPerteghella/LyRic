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

    // Call the process or workflow from helpers.nf
    fastqStats()
}