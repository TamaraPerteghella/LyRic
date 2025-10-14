process makeStarIndex{
    tag "Build Genome index"

    input:
    tuple val(file_name), val(genome)

    output:
    path "SA/"
    
    script:
    """
    #mkdir -p "${params.genomes_folder}/STARshort_indices/${params.genome}/SA/"
    #mkdir SA/
    STAR --runMode genomeGenerate --runThreadN 3 --genomeDir "SA/" --genomeFastaFiles "${params.genomes_folder}/${params.genome}.fa.gz"
    """
}

//STILL AWAITING TO BE COMPLETE
process hiSeqReadMapping{
    input:
    val(file_name)

    output:
    tuple path("hiSeq_${file_name}.bam"), path("hiSeq_${file_name}.bam.bai")
    
    script:
    """
    #mkdir -p "${params.mappingsdir}/shortReadMappings/" "${params.mappingsdir}/STAR/${file_name}"
    
    echoerr "Mapping"

    STAR --runThreadN ${task.cpus} --readFilesIn ${params.shortreadsdir}/${file_name}_1.fastq.gz ${params.shortreadsdir}/${file_name}_2.fastq.gz \
        --genomeDir "${params.genomes_folder}/STARshort_indices/${params.genome}/SA/" \
        --readFilesCommand zcat --sjdbOverhang 124 \
        --outFileNamePrefix "output/mappings/STAR/`basename {output}`/" \
        --outStd SAM --genomeLoad NoSharedMemory \
        --outSAMunmapped Within --outFilterType BySJout \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
        --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS nM XS \
    | samtools view -b -u -S - | samtools sort -T ${params.TMPDIR} -@ 2 -m 15000000000 - > hiSeq_${file_name}.bam

    sleep 200s
    samtools index hiSeq_${file_name}.bam
    echoerr "Mapping done"
   """
}

workflow srMapping {
    take:
    input_samples

    main:
    if ( "${params.short_read}".toBoolean() ) {
        println("Proceeding with short'read mapping with STAR")

        index = makeStarIndex( input_samples )



    } else {
        println("Short read mapping - SKIPPED")
        index = ""
    }

    emit:
    index

}