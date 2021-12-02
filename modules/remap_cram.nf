
//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//    log.info "${params.ref_dir}"

process remap_cram {
    //memory '12G'
    tag "$sample"
    cpus 4
    disk '20 GB'
    //time '100m'
    //queue 'normal'
    clusterOptions = { "-M 18000 -R \"select[mem>=18000] rusage[mem=18000]\" -n4 -R \"span[hosts=1]\" -R \"select[model==Intel_Platinum]\"" }
    //container  = 'file:///software/hgi/containers/samtools_sambamba.sif'
    //containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { (task.attempt <= maxRetries)  ? 'retry' : 'ignore' }
    //publishDir "${params.outdir}/cram_index/", mode: 'symlink', overwrite: true, pattern: "${cram_file}.crai"
    maxRetries 3

    when:
    params.remap
    //params.run
     
    input:
    tuple val(study_id), val(sample), path(cram_file), path(cram_file_index)
    //path cram_file

    output:
    tuple val(study_id), val(sample), path("${sample}.remapped.cram"), emit: remapped_sample_cram
    //path "${cram_file}.sorted"
    //tuple file("${cram_file}.sorted"), emit: indexes

    script:
""" 
singularity run -B ${params.ref_dir}:${params.ref_dir} -B \$(pwd):/home --pwd "/home" --no-home /software/hgi/containers/oqfe_remap.sif -1 /home/${cram_file} --sample ${sample} --cram-reference-fasta ${params.ref_dir}/hs38DH.fa -c -j 4
"""
}

