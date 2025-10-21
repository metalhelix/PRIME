include {run_base2fastq; namechange_and_copy_fastqs; run_cellranger; combined_report_generate} from "../modules/Aviti_core_10X.nf"
include {SAMPLESHEET} from "../workflows/SampleSheetAviti.nf"
include {DRIVERCSV} from "./DriverG4.nf"
include {driver_make; bowtie2_fastqc; bam_stats; flagstat; copy_fastqs; multiQC; sample_report_make; multiQC_report_provide} from "../modules/G4_core_RNAseq.nf"
include {bowtie2_fastqc_Aviti; driver_make_Aviti} from "../modules/Aviti_core_RNAseq.nf"
include {cellbender; seurat_report; shiny_app_make; biotools_append} from "../modules/CellBender_SeuratInit.nf"

workflow AVITI_10X {
    take:
    lims_info_csv

    main:
    meta = SAMPLESHEET(lims_info_csv)
    // Determine orgin of the fastq directory
    if (!params.fastq_dir) {
        run_base2fastq(lims_info_csv)
        fastqs = run_base2fastq.out
    } else {
        fastqs = params.fastq_dir
    }
    // 
    namechange_and_copy_fastqs(meta, fastqs)
    run_cellranger(namechange_and_copy_fastqs.out, lims_info_csv)
    cellbender(run_cellranger.out)
    seurat_report(cellbender.out)

    def all_seurat_report_out = seurat_report.out
                                            .collect()
                                            .map{it[0,1,2]}

    combined_report_generate(all_seurat_report_out)
    shiny_app_make(all_seurat_report_out, combined_report_generate.out)
    biotools_append(all_seurat_report_out, combined_report_generate.out)

    emit:
    combined_report_generate.out
}

workflow AVITI_RNASEQ { 
    take:
    lims_info_csv

    main:
    meta = SAMPLESHEET(lims_info_csv)
   // Determine orgin of the fastq directory
    if (!params.fastq_dir) {
        run_base2fastq(lims_info_csv)
        fastqs = run_base2fastq.out
    } else {
        fastqs = params.fastq_dir
    }
    //  
    namechange_and_copy_fastqs(meta, fastqs)
    def all_nameChange_outs = namechange_and_copy_fastqs.out
                                                        .collect()
                                                        .map {it[0]}
    driver_csv = driver_make_Aviti(lims_info_csv, fastqs, all_nameChange_outs) 
    driver = DRIVERCSV(driver_csv)
    bowtie2_fastqc(driver)
    bam_stats(bowtie2_fastqc.out)
    flagstat(bowtie2_fastqc.out)
    
    def all_bowtieOuts = bowtie2_fastqc.out
                                    .collect()
                                    .map {it[0]}

    multiQC(all_bowtieOuts, driver_csv)
    sample_report_make(driver_csv, multiQC.out)
    multiQC_report_provide(sample_report_make.out)

    emit:
    multiQC_report_provide.out

}
