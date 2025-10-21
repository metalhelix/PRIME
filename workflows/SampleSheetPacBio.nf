def create_metadata_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.libID = row.LibID
    meta.nanalysis_path = row.resultPaths
    meta.molngID = row.Molng
    meta.bam_path = row.BamPath
    meta.barcode = row.Barcode
    meta.qcReport = row.QcReport
    meta.fcpath = row.FcPath

    return meta
    }

process samplesheet_make_pacBio {

    label 'small_mem'

    input:
    val pacbio_xml

    output:
    path "pacbio_metadata.csv"

    script:
    """
    ml python/3.11.5 

    pacbio_metaData_reader.py --xml ${pacbio_xml}
    """
}

workflow SAMPLESHEET_PACBIO {
    take:
    pacbio_xml

    main:
    samplesheet_make_pacBio(pacbio_xml)
    samplesheet_make_pacBio.out.splitCsv ( header:true, sep:',' )
        .map { create_metadata_channel(it) }
        .set { data_meta }
    
    emit:
    data_meta
}