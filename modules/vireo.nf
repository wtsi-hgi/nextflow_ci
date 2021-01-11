process vireo {
    tag "${samplename}"
    publishDir "${params.outdir}/vireo/", mode: 'copy', overwrite: true
    
    input:
    tuple val(samplename), path(cell_data), val(n_pooled)
    
    output:
    tuple val(samplename), path("vireo_${samplename}"), emit: vireo_output_dir

    script:
    """
umask 2 # make files group_writable

vireo -c $cell_data -N $n_pooled -o vireo_${samplename}
    """
}
