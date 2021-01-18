process plot_donor_ncells {
    tag ""
    
    publishDir "${params.outdir}/", mode: 'copy', pattern: "*.pdf", overwrite: true
    
    when: 
    params.plot_donor_ncells.run

    input: 
    path(sample_donor_summary_tsv)

    output: 
    path("*.pdf"), emit: pdf

    script:
    """
    python $workflow.projectDir/../bin/plot_donor_ncells.py ${sample_donor_summary_tsv}
    """
}
