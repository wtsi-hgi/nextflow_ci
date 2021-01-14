process split_donor_h5ad {
    tag "${sample} ${run_id} ${study_id}"
    
    when: 
    params.run_imeta_study_cellranger

    input: 
    tuple val(study_id), val(sample), val(run_id)

    output: 
    tuple val(study_id), val(sample), val(run_id), env(CELLRANGER_IRODS_OBJECT), env(WORK_DIR), emit: study_id_sample_cellranger_object
    env(WORK_DIR), emit: work_dir_to_remove

    script:
    """
    python $workflow.projectDir/../bin/split_h5ad_per_donor.py ${sample} ${run_id}
    """
}