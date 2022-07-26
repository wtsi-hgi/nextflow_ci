process tximport {
    tag "$params.ensembl_lib"
    
    publishDir "${params.outdir}/tximport", mode: "${params.copy_mode}", overwrite: true

    when:
    params.salmon_aligner.salmon_dowstream_tasks.tximport_task.run

    input:
    file (quant_sf_files)  // from collect()

    output:
    file("fofn_quantfiles.txt")
    file("txi_gene_counts.csv")
    file("txi_transcript_counts.csv")
    file("txi_lengthScaledTPM_gene_counts.csv")
    file("tximport.rdata")
    //file "${samplename}.quant.sf" // into ch_salmon_trans
    //file "${samplename}.quant.genes.sf" //into ch_salmon_genes
    // file "my_outs/${samplename}" optional true // into ch_alignment_logs_salmon

    script:
    """
    ls . | grep .quant.sf\$ > fofn_quantfiles.txt

    /software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript $workflow.projectDir/../bin/tximport.R \"${params.salmon_aligner.salmon_downstream_tasks.tximport_task.ensembl_lib}\" fofn_quantfiles.txt 
    """
}
