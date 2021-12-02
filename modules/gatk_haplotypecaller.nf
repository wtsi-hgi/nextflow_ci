//params.run = true

//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//params.bed_dir = "/lustre/scratch118/humgen/resources/exome/Homo_sapiens/"

process gatk_haplotypecaller {
    //memory '8G'
    tag "$sample"
    cpus 2
    //disk '20 GB'
    //time '100m'
    //queue 'normal'
    clusterOptions = { "-M 8000 -R \"select[mem>=8000] rusage[mem=8000]\" -R \"select[model==Intel_Platinum || model==Intel_Platinuma || model==Intel_Platinumb || model==Intel_Platinumc || model==Intel_Platinumd || model==Intel_Platinume || model==Intel_Platinumf || model==Intel_Platinumg]\"" }
    container  = 'file:///software/hgi/containers/gatk_4.2.2.0.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind ${params.bed_dir}:/bed_files --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { (task.attempt <= maxRetries)  ? 'retry' : 'ignore' }
    publishDir "${params.haplotypecaller_output_dir}", mode: 'copy', overwrite: true, pattern: "*gatk.g.vcf.gz*"
    maxRetries 3

    when:
    params.run_haplotypecaller
     
    input:
    tuple val(study_id), val(sample), path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)

    output:
    tuple path("${cram_file_sorted_dups_coord}.gatk.g.vcf.gz"), path("${cram_file_sorted_dups_coord}.gatk.g.vcf.gz.tbi")
    //tuple file("${cram_file}.sorted"), emit: indexes


    script:
""" 
/gatk/gatk --java-options "-Xms6g -Xmx6g  -XX:+UseSerialGC" HaplotypeCaller -I ${cram_file_sorted_dups_coord} -O ${cram_file_sorted_dups_coord}.gatk.g.vcf.gz -R /ref/hs38DH.fa -L /bed_files/Twist/Twist_Human_Core_Exome_BI-CTR_padded_merged.interval_list -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 -ERC GVCF -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation
"""
}

