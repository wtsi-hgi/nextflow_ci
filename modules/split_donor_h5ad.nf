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
    bash $workflow.projectDir/../bin/imeta_study_cellranger.sh ${sample} ${run_id}
    if [ -f cellranger.object.txt ] 
    then 
        echo file cellranger.object.txt found
        CELLRANGER_IRODS_OBJECT=\$(cat cellranger.object.txt)
        WORK_DIR=dont_remove
        rm cellranger.object.txt
    else
        echo not found file cellranger.object.txt
        CELLRANGER_IRODS_OBJECT=cellranger_irods_not_found
        WORK_DIR=\$PWD
    fi
    """
}
