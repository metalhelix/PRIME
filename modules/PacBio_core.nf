process copy_all_files {
    label 'small_mem'

    input:
    val meta

    output:
    val "${meta.nanalysis_path}"

    script:
    """
    mkdir -p ${meta.nanalysis_path}

    rsync ${meta.bam_path} ${meta.nanalysis_path}/${meta.libID}.bam
    rsync ${meta.qcReport} ${meta.nanalysis_path}/SMRT_QCReport.pdf

    rsync -a --ignore-existing ${meta.fcpath}/fail_reads ${meta.nanalysis_path}
    rsync -a --ignore-existing ${meta.fcpath}/hifi_reads/*.unassigned.bam ${meta.nanalysis_path}/unassigned.bam
    """
}

process pacbio_samplereport {
    label 'small_mem'

    input:
    val nanalysis_path
    path samplesheet_pacbio

    publishDir "${nanalysis_path}" , mode: 'copy'

    output:
    path "Sample_Report.csv"
    val "https://webfs/${nanalysis_path}Sample_Report.csv ", emit: sample_report_link
    val "https://webfs/${nanalysis_path}SMRT_QCReport.pdf ", emit: qc_report_link

    script:
    """
    ml python/3.11.5

    Sample_Report_PacBio.py --input ${samplesheet_pacbio}
    """
}