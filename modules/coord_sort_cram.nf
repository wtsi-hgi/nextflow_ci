//params.run = true

//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process coord_sort_cram {
    //memory '12G'
    tag "$cram_file"
    //cpus 1
    disk '20 GB'
    //time '100m'
    //queue 'normal'
    clusterOptions = { "-M 18000 -R \"select[mem>=18000] rusage[mem=18000]\" -R \"select[model==Intel_Platinum]\"" }
    container  = 'file:///software/hgi/containers/sambamba_0.6.4.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { (task.attempt <= maxRetries)  ? 'retry' : 'ignore' }
    //publishDir "${params.outdir}/cram_index/", mode: 'symlink', overwrite: true, pattern: "${cram_file}.crai"
    maxRetries 3

    when:
    params.run_coord_sort_cram
     
    input:
    tuple val(study_id), val(sample), path(cram_file_sorted_dups)
    //path cram_file_sorted_dups

    output:
    tuple val(study_id), val(sample), path("${cram_file_sorted_dups}.coord"), path("${cram_file_sorted_dups}.coord.bai"), emit: markdup_sample_cram_crai
    //tuple path("${cram_file_sorted_dups}.coord"), path("${cram_file_sorted_dups}.coord.bai")
    //tuple file("${cram_file}.sorted"), emit: indexes

    script:
""" 
sambamba sort -p -m 7GB --tmpdir ./tmp ${cram_file_sorted_dups} -o ${cram_file_sorted_dups}.coord && sambamba index ${cram_file_sorted_dups}.coord
"""
}

