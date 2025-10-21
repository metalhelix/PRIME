include {SAMPLESHEET_PACBIO; samplesheet_make_pacBio} from "./SampleSheetPacBio.nf"
include {copy_all_files; pacbio_samplereport} from "../modules/PacBio_core.nf"

workflow PACBIO {
    take:
    pacbio_xml

    main:
    SAMPLESHEET_PACBIO(pacbio_xml)

    meta = SAMPLESHEET_PACBIO.out
    samplesheet_csv = samplesheet_make_pacBio(pacbio_xml)

    copy_all_files(meta)

    nanalysis_path = copy_all_files.out.collect().map {it[0]}

    pacbio_samplereport(nanalysis_path, samplesheet_csv)

    emit:
    sample_report_link = pacbio_samplereport.out.sample_report_link
    qc_report_link = pacbio_samplereport.out.qc_report_link
}