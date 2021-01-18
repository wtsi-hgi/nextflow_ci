process plot_donor_ncells {
    tag "${sample_donor_summary_tsv}"
    
    publishDir "${params.outdir}/", mode: "${params.copy_mode}", pattern: "*.pdf", overwrite: true
    
    when: 
    params.plot_donor_ncells.run

    input: 
    path(sample_donor_summary_tsv)

    output: 
    path("*.pdf"), emit: pdf

    script:
    """
python $workflow.projectDir/../bin/plot_donor_ncells.py \\
  --sample_donor_summary_tsv ${sample_donor_summary_tsv} \\
  --plotnine_dpi ${params.plotnine_dpi.plotnine_dpi}
    """
}
