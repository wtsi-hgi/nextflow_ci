//params.run = true

//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//params.bed_dir = "/lustre/scratch118/humgen/resources/exome/Homo_sapiens/"

process deepvariant {
    //memory '8G'
    tag "$sample"
    cpus 2
    //disk '20 GB'
    //time '100m'
    //queue 'normal'
    clusterOptions = { "-M 8000 -R \"select[mem>=8000] rusage[mem=8000]\" -R \"select[model==Intel_Platinum]\"" }
    container  = 'file:///software/hgi/containers/deepvariant_0.10_UKB.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind ${params.bed_dir}:/bed_files --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { (task.attempt <= maxRetries)  ? 'retry' : 'ignore' }
    publishDir "${params.deepvariant_output_dir}/gvcf", mode: 'copy', overwrite: true, pattern: "*dv.g.vcf.gz*"
    publishDir "${params.deepvariant_output_dir}/vcf", mode: 'copy', overwrite: true, pattern: "*dv.vcf.gz*"
    maxRetries 3

    when:
    params.run_deepvariant
     
    input:
    tuple val(study_id), val(sample), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)

    output:
    tuple path("${cram_file_sorted_dups_coord}.dv.vcf.gz"), path("${cram_file_sorted_dups_coord}.dv.vcf.gz.tbi"), path("${cram_file_sorted_dups_coord}.dv.g.vcf.gz"), path("${cram_file_sorted_dups_coord}.dv.g.vcf.gz.tbi")
    //tuple file("${cram_file}.sorted"), emit: indexes


    script:
""" 
/opt/deepvariant/bin/run_deepvariant --model_type=WES --customized_model=/opt/models/ukb_wes/model.ckpt-22236 --ref=/ref/hs38DH.fa --reads=${cram_file_sorted_dups_coord} --output_vcf=${cram_file_sorted_dups_coord}.dv.vcf.gz --output_gvcf=${cram_file_sorted_dups_coord}.dv.g.vcf.gz --intermediate_results_dir ./tmp --num_shards=1 --regions=/bed_files/Twist/Twist_Human_Core_Exome_BI-CTR_padded_merged.bed
"""
}

