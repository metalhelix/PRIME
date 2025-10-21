process cellbender {
    label 'gpu'
    // conda "${projectDir}/env/cellbender.yml"
    conda "/home/compbio_svc/miniconda3/envs/cellbender"

    input:
    tuple val(secondary_path), val(orderType), path(cellranger_outs), val(output_lib_folder), val(meta)

    output:
    tuple val(secondary_path), val(orderType), val(output_lib_folder), val(meta), path("cellbender_output*")

    publishDir "${output_lib_folder}/cellbender_outs", mode: 'copy'

    script:
    """
    mkdir -p ${output_lib_folder}/cellbender_outs

    cellbender remove-background \
    --input  ${cellranger_outs}/raw_feature_bc_matrix.h5 \
    --output cellbender_output.h5 \
    --learning-rate 0.0001 \
    --cuda

    ptrepack --complevel 5 cellbender_output_filtered.h5:/matrix cellbender_output_filtered_seurat.h5:/matrix
    """
}

process seurat_report {
    label 'big_mem'
    // conda "${projectDir}/env/SC_PRIME.yml"
    conda "/home/compbio_svc/miniconda3/envs/SC_PRIME"

    input:
    tuple val(secondary_path), val(orderType), val(output_lib_folder), val(meta), val(cellbender_output_h5file)

    output:
    tuple val(secondary_path), val(orderType), val(output_lib_folder)

    script:

    if (!params.annotation) {
        annotation = meta.annotation
    } else {
        annotation = params.annotation
    }

    if (!params.genome){
        index_genome = "${meta.species}/${meta.refGenome}"
        genome_ver = meta.refGenome
    } else {
        index_genome = params.genome
        genome_ver = index_genome.split("/")[1]
    }

    """
    Seurat_ReportGen.py --cellbender_folder ${output_lib_folder}/cellbender_outs --rmd ${projectDir}/assets/Seurat_Report_gen.rmd --mt_file ${params.indexDir}/${index_genome}/annotation/${annotation}/extras/${genome_ver}.${annotation}.mtGeneList.txt
    """
}

process shiny_app_make {
    label 'lil_mem'

    input:
    tuple val(secondary_path), val(orderType), val(output_lib_folder)
    val report_link
    
    output:
    tuple val(secondary_path), val(orderType), val(output_lib_folder)

    script:
    """
    Shiny_app_deploy.py --second_folder ${secondary_path} --wwwPath ${projectDir}/assets/SC_shiny/www
    """
}

process biotools_append {
    label 'small_mem'

    input:
    tuple val(secondary_path), val(orderType), path(cellranger_outs) 
    val report_link

    output:
    val "Uploaded"

    publishDir "${params.fcpath}", mode: 'copy'

    script:
    """
    mongodb_append.py -s ${secondary_path}
    """
}
