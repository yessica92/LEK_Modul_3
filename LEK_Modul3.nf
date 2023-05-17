nextflow.enable.dsl = 2

params.accession = "null"
params.outdir = "SRA_data"

process prefetch {
  storeDir "${params.outdir}"
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

  input:
    val accession

  output:
    path "sra/${accession}.sra", emit: srafile 

  script:
    """
    prefetch $accession 
    """
}

process convertfastq {
  storeDir "${params.outdir}/sra/${accession}_fastq/raw_fastq"

  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

  input: 
    val accession
    path srafile

  output:
     path '*.fastq', emit: fastqfiles

  script:
  """
  fastq-dump --split-files ${srafile}
  """
}

process fastqc {
    publishDir "${params.outdir}/sra/${params.accession}_fastq/", mode: "copy", overwrite: true
    container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"

    input:
      path fastqfiles

    output:
      path "fastqc_results/*.html", emit: results
      path "fastqc_results/*.zip", emit: fastqczip
      

    script:
      """
      mkdir fastqc_results
      fastqc ${fastqfiles} --outdir fastqc_results
      """
}

process fastp {
  publishDir "${params.outdir}/sra/${accession}_fastq/", mode: 'copy', overwrite: true
  
  container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"

  input:
    path fastqfiles
    val accession

  output:
    path "fastp_fastq/*.fastq", emit: fastqfiles
    path "fastp_report", emit: fastpreport

  script:
      if(fastqfiles instanceof List) {
      """
      mkdir fastp_fastq
      mkdir fastp_report
      fastp -i ${fastqfiles[0]} -I ${fastqfiles[1]} -o fastp_fastq/${fastqfiles[0].getSimpleName()}_fastp.fastq -O fastp_fastq/${fastqfiles[1].getSimpleName()}_fastp.fastq -h fastp_report/fastp.html -j fastp_report/fastp.json
      """
    } else {
      """
      mkdir fastp_fastq
      mkdir fastp_report
      fastp -i ${fastqfiles} -o fastp_fastq/${fastqfiles.getSimpleName()}_fastp.fastq -h fastp_report/fastp.html -j fastp_report/fastp.json
      """
    }
}


workflow {
    srafile = prefetch(Channel.from(params.accession))
    converted = convertfastq(Channel.from(params.accession), srafile)
    converted_flat = conveerted.flatten()
    reports_channel = Channel.empty()
    fastqc_outputchannel = fastqc(converted_flat)
    reports_channel = reports_channel.concat(fastqc_outputchannel)
    reports_channel_collected = reports_channel.collect()
    reports_channel_collected.view()
}