nextflow.preview.dsl=2

params.runtag = 'deepvariant'
//params.index_crams = false
//params.run_deepvariant = true
//params.outdir = "${baseDir}/../../outputs"
//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//params.bed_dir = "/lustre/scratch118/humgen/resources/exome/Homo_sapiens/"

include { sort_cram } from '../modules/sort_cram.nf' //params(run: true, outdir: params.outdir)
include { markDuplicates } from '../modules/markDuplicates.nf' //params(run: true, outdir: params.outdir)
include { coord_sort_cram } from '../modules/coord_sort_cram.nf' //params(run: true, outdir: params.outdir)
include { deepvariant } from '../modules/deepvariant.nf' //params(run: true, outdir: params.outdir)

workflow {

    log.info "${params}"

    if (params.run_deepvariant) {
	

	ch_cram_file = Channel
		.fromPath(params.cram_fofn)
		.splitText()
		.take(1)
		//.view()
     main:
        
        log.info "${params.ref_dir}"
        
	sort_cram(ch_cram_file)
	markDuplicates(sort_cram.out)
	coord_sort_cram(markDuplicates.out)
	deepvariant(coord_sort_cram.out)
\\     emit:
\\        my_data = deepvariant.out
        
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
