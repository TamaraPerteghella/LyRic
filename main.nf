#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
nextflow.preview.output = true

// Import everything from helpers.nf
include { fastqStats } from './processes/fastqStats.nf'
include { lrMapping } from './processes/lrMapping.nf'
include { srMapping } from './processes/srMapping.nf'

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
                def (filePath, filename, genome, tech) = line.split('\t')[0..3]
                tuple(filename, file(filePath), genome, tech)
            }
        fastq_ch.view()
    }
    else {
        println("Expecting files in: ${params.datadir}/")
        fastq_ch = Channel.fromPath("${params.datadir}/*.fastq.gz")
            .map { file ->
                def base = file.name.replaceAll(/\.fastq\.gz$/, '')
                tuple(base, file, "${params.genome}", "${params.tech}")
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
    fastqStats(fastq_ch.map{ file_name, fastq, _genome, _tech -> tuple(file_name, fastq) })

    statfiles = fastqStats.out.timestamp_fastq.concat(fastqStats.out.aggregate_summary)
    qcfiles = fastqStats.out.qc_ch
    plotfiles = fastqStats.out.read_length_plot
    tmpstatsfile = fastqStats.out.readlength_ch.map { t -> t[1] }.concat(fastqStats.out.readlength_summary_ch)

    lrMapping(fastq_ch)

    long_reads = lrMapping.out.mappings.concat(lrMapping.out.indexes)
    bwigs = lrMapping.out.bigwigs

    lr_bamqc = lrMapping.out.bamqc_ch
    lr_qc = lrMapping.out.dupl

    statfiles = statfiles.concat(lrMapping.out.agg_stats).concat(lrMapping.out.matrix).concat(lrMapping.out.allbasic).concat(lrMapping.out.allspikes)
    plotfiles = plotfiles.concat(lrMapping.out.plots).concat(lrMapping.out.density).concat(lrMapping.out.heatmap).concat(lrMapping.out.plot_stats).concat(lrMapping.out.plot_spikeins)

    exonic_bw = lrMapping.out.exonic_bigwigs
    lr_bed = lrMapping.out.beds
    lr_gffs = lrMapping.out.gffs

    biotype_class = lrMapping.out.biotype_class
    statfiles = statfiles.concat(lrMapping.out.biotype_stats)
    plotfiles = plotfiles.concat(lrMapping.out.plot_biotype_stats)

    srMapping(fastq_ch.map{ file_name, _fastq, genome, _tech -> tuple(file_name, genome) })

    publish:
    qc = qcfiles
    stats = statfiles
    stats_tmp = tmpstatsfile
    plots = plotfiles
    long_reads_mappings = long_reads
    long_reads_bwig = bwigs
    long_reads_bamqc = lr_bamqc
    long_reads_exonic_bwig = exonic_bw
    long_reads_qc = lr_qc
    long_reads_bed = lr_bed
    long_reads_gff = lr_gffs
    btp = biotype_class
}

output {
    qc {
        path { _fastqc -> "${params.qcdir}/" }
    }
    plots {
        path { _plt -> "${params.plotdir}/" }
    }
    stats {
        path { _stat -> "${params.statsdir}/" }
    }
    stats_tmp {
        path { _stat -> "${params.statsdir}/tmp/" }
    }
    long_reads_mappings {
        path { _lrm -> "${params.longmappingsdir}" }
    }
    long_reads_bwig {
        path { _lrb -> "${params.longmappingsdir}/bigwigs/" }
    }
    long_reads_bamqc {
        path { _lrq -> "${params.longmappingsdir}/qc/bamqc/" }
    }
    long_reads_exonic_bwig {
        path { _lre -> "${params.longmappingsdir}/exonic_bigwigs/" }
    }
    long_reads_qc {
        path { _lrq -> "${params.longmappingsdir}/qc/" }
    }
    long_reads_bed {
        path { _lrb -> "${params.longmappingsdir}/readBamToBed/" }
    }
    long_reads_gff {
        path { _lrg -> "${params.longmappingsdir}/readBedToGff/" }
    }
    btp {
        path { _btp -> "${params.longmappingsdir}/reads2biotypes/" }
    }
}
