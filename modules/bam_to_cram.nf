//params.run = true

//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process bam_to_cram {
    //memory '12G'
    tag "$sample"
    //cpus 1
    disk '20 GB'
    //time '100m'
    //queue 'normal'
    clusterOptions = { "-M 18000 -R \"select[mem>=18000] rusage[mem=18000]\" -R \"select[model==Intel_Platinuma || model==Intel_Platinumb || model==Intel_Platinumc || model==Intel_Platinumd || model==Intel_Platinume || model==Intel_Platinumf || model==Intel_Platinumg]\"" }
    container  = 'file:///software/hgi/containers/samtools-1.10.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { (task.attempt <= maxRetries)  ? 'retry' : 'ignore' }
    //publishDir "${params.outdir}/cram_index/", mode: 'symlink', overwrite: true, pattern: "${cram_file}.crai"
    publishDir "${params.cram_output_dir}", mode: "${params.copy_mode}"
    maxRetries 3

    when:
    params.run_coord_sort_cram
     
    input:
    tuple val(study_id), val(sample), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)
    //path cram_file_sorted_dups

    output:
    tuple val(study_id), val(sample), path("${sample}.sorted.dups.coord.cram"), path("${sample}.sorted.dups.coord.cram.crai"), emit: processed_sample_cram_crai
    //tuple path("${cram_file_sorted_dups}.coord"), path("${cram_file_sorted_dups}.coord.bai")
    //tuple file("${cram_file}.sorted"), emit: indexes

    script:
""" 
samtools view ${cram_file_sorted_dups_coord} -o ${sample}.sorted.dups.coord.cram -O CRAM -T ${params.ref_dir}/hs38DH.fa && samtools index ${sample}.sorted.dups.coord.cram
"""
}

