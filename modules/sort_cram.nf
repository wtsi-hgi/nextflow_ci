
//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//    log.info "${params.ref_dir}"

process sort_cram {
    memory '8G'
    tag "$cram_file"
    //cpus 1
    disk '20 GB'
    //time '100m'
    //queue 'normal'
    container  = 'file:///software/hgi/containers/sambamba_0.6.4.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    //publishDir "${params.outdir}/cram_index/", mode: 'symlink', overwrite: true, pattern: "${cram_file}.crai"
    maxRetries 3

    when:
    params.run_sort_cram
    //params.run
     
    input:
    path cram_file

    output:
    path "${cram_file}.sorted"
    //tuple file("${cram_file}.sorted"), emit: indexes

    script:
""" 
/software/sciops/pkgg/samtools/1.10+42_g3c4b380+irods_4.2.7/bin/samtools view -O BAM ${cram_file} | sambamba sort -p -n --tmpdir /tmp /dev/stdin -o ${cram_file}.sorted
"""
}

