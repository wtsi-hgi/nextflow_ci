nextflow.enable.dsl=2

// All inputs are read from Nextflow config file "inputs.nf",
//  which is located in upstream Gitlab "nextflow_ci" repo (at same branch name).
// Meaning that if you wish to run pipeline with different parameters,
// you have to edit+commit+push that "inputs.nf" file, then rerun the pipeline.

// import modules that depend on input mode:
include { test_lustre_access } from '../modules/test_lustre_access.nf'
include { cellsnp } from '../modules/cellsnp.nf'
include { vireo } from '../modules/vireo.nf'

workflow {

    if (params.test_lustre_access.run_test) {
	test_lustre_access(Channel.from(params.test_lustre_access.lustre_dir))
	// work_dir_to_remove = imeta_study.out.work_dir_to_remove
    }
    
    if (params.cellsnp.run) {

	Channel.fromPath(params.cellsnp.cellranger_input.lustre_filepath10x_tsv)
            .splitCsv(header: true, sep: '\t')
	    .map{row->tuple(row.experiment_id, row.data_path_10x_format)}
	    .set{ch_experiment_path10x}

	if (params.cellsnp.cellranger_input.replace_lustre_path) {
	    ch_experiment_path10x
		.map{experiment, path10x ->
		tuple(experiment,
		      path10x.replaceFirst(/${params.cellsnp.cellranger_input.replace_path_from}/,
					   params.cellsnp.cellranger_input.replace_path_to))}
		.set{ch_experiment_path10x_tocellsnp}
	} else {
	    ch_experiment_path10x
		.set{ch_experiment_path10x_tocellsnp}
	    
	}

	cellsnp(, Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())
	// work_dir_to_remove = imeta_study.out.work_dir_to_remove
    }
    
    if (params.vireo.run) {
    }
}

workflow.onError {
    log.info "Pipeline execution stopped with the following message: ${workflow.errorMessage}" }

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    
    if (params.on_complete_uncache_irods_search) {
	log.info "You have selected \"on_complete_uncache_irods_search = true\"; will therefore attempt to remove Irods work dirs to forcefully uncache them even if successful."
	if (! file("${params.outdir}/irods_work_dirs_to_remove.csv").isEmpty()) {
	    log.info "file ${params.outdir}/irods_work_dirs_to_remove.csv exists and not empty ..."
	    file("${params.outdir}/irods_work_dirs_to_remove.csv")
		.eachLine {  work_dir ->
		if (file(work_dir).isDirectory()) {
		    log.info "removing work dir $work_dir ..."
		    file(work_dir).deleteDir()   
		} } } }
    
    if (params.on_complete_remove_workdir_failed_tasks) {
	log.info "You have selected \"on_complete_remove_workdir_failed_tasks = true\"; will therefore remove work dirs of all tasks that failed (.exitcode file not 0)."
	// work dir and other paths are hardcoded here ... :
	def proc = "bash ./nextflow_ci/bin/del_work_dirs_failed.sh ${workDir}".execute()
	def b = new StringBuffer()
	proc.consumeProcessErrorStream(b)
	log.info proc.text
	log.info b.toString() }
}


	//        .filter { it[2] =~ /.cram$/ } // Need to check for bam too?
	//        .unique()
