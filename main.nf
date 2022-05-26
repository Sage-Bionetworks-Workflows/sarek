#!/usr/bin/env nextflow

tsvFile = file(params.input)
inputSample = extractFastq(tsvFile)

(genderMap, statusMap, inputSample) = extractInfos(inputSample)

ch_fasta = Channel.value(file(params.fasta))
ch_fai   = Channel.value(file(params.fasta_fai))
ch_bwa   = Channel.value(file(params.bwa))

inputPairReads = inputSample

(inputPairReads, inputPairReadsFastQC) = inputPairReads.into(2)

process FastQCFQ {
    label 'FastQC'
    label 'cpus_2'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsFastQC

    output:
        file("*.{html,zip}") into fastQCFQReport

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}

process MapReads {
    label 'cpus_max'

    tag "${idPatient}-${idRun}"

    input:
        set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from inputPairReads
        file(bwaIndex) from ch_bwa
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bamMapped
        set idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam") into bamMappedBamQC

    script:
    status = statusMap[idPatient, idSample]
    """
    bwa mem -K 100000000 -t ${task.cpus} -M ${fasta} ${inputFile1} ${inputFile2} | \
    samtools sort --threads ${task.cpus} -m 2G - > ${idSample}_${idRun}.bam
    """
}

def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    // The next line forces roundtrip as a list to ensure
    // that genderMap and statusMap are fully populated
    // before being returning and used in processes
    channel = Channel.fromList(channel.toList().val)
    [genderMap, statusMap, channel]
}

def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = row[2].toInteger()
            def idSample   = row[3]
            def idRun      = row[4]
            def file1      = returnFile(row[5])
            def file2      = returnFile(row[6])
            [idPatient, gender, status, idSample, idRun, file1, file2]
        }
}

def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}
