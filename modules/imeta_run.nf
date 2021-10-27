process imeta_run {
    tag "${id_run}"
    publishDir "${params.outdir}/imeta_run/run_id_${id_run}/", mode: 'copy', pattern: "samples.tsv", overwrite: true
    publishDir "${params.outdir}/imeta_run/run_id_${id_run}/", mode: 'copy', pattern: "samples_noduplicates.tsv", overwrite: true
    publishDir "${params.outdir}/", mode: 'copy', pattern: "samples.tsv", saveAs: { filename -> "${id_run}.$filename" }, overwrite: true
    publishDir "${params.outdir}/", mode: 'copy', pattern: "samples_noduplicates.tsv", saveAs: { filename -> "${id_run}.$filename" }, overwrite: true

    when: 
    params.run_id_mode.run_imeta_run

    input: 
    val(id_run)

    output: 
    tuple val(id_run), path('samples.tsv'), emit: irods_samples_tsv
    tuple val(id_run), path('samples_noduplicates.tsv'), emit: samples_noduplicates_tsv
    tuple val(id_run), path('guides_library.tsv'), emit: guides_library_tsv
    env(WORK_DIR), emit: work_dir_to_remove

    script:
    """
    bash $workflow.projectDir/../bin/imeta_run.sh ${id_run} ${params.guide_library}
    awk '!a[\$1]++' samples.tsv > samples_noduplicates.tsv 

    # Save work dir so that it can be removed onComplete of workflow, 
    # to ensure that this task Irods search is re-run on each run NF run, 
    # in case new sequencing samples are ready: 
    WORK_DIR=\$PWD
    """
}
// awk removes duplicates as one sanger sample can have several run_id
