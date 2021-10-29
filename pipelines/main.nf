nextflow.enable.dsl=2

// All inputs are read from Nextflow config file "inputs.nf",
//  which is located in upstream Gitlab "nextflow_ci" repo (at same branch name).
// Meaning that if you wish to run pipeline with different parameters,
// you have to edit+commit+push that "inputs.nf" file, then rerun the pipeline.

params.guide_libraries = "${baseDir}/../../inputs/guide_libraries/sgRNA_Minlib_CTLA4.guide_library.csv"

Channel.fromPath(params.guide_libraries)
    .set{ch_library_files}

params.samplename_library = "${baseDir}/../../inputs/sample_guides.csv"
Channel.fromPath(params.samplename_library)
    .splitCsv(header: true)
    .map { row -> tuple("${row.samplename}", "${row.library}", "${row.includeG}") }
    .set{ch_samplename_library}

// import modules that depend on input mode:
include { imeta_study } from '../modules/imeta_study.nf'
include { imeta_run } from '../modules/imeta_run.nf'
include { imeta_samples_csv } from '../modules/imeta_samples_csv.nf'
include { gsheet_to_csv } from '../modules/gsheet_to_csv.nf'
// module specific to google_spreadsheet input mode:
include { join_gsheet_metadata } from '../modules/join_gsheet_metadata.nf'
include { iget_study_cram } from '../modules/iget_study_cram.nf'
include { crams_to_fastq } from '../modules/crams_to_fastq.nf'
include { fastx_trimmer } from '../modules/crispr/trim_fastq.nf' params(run: true, outdir: params.outdir)
include { merge_fastq_batches } from '../modules/crispr/merge_fastq_batches.nf' params(run:true, outdir: params.outdir)
include { count_crispr_reads } from '../modules/crispr/count_crispr_reads.nf' params(run: true, outdir: params.outdir,
                                                         read2: params.read2)
include { collate_crispr_counts } from '../modules/crispr/collate_crispr_counts.nf' params(run: true, outdir: params.outdir)
include { fastqc } from '../modules/crispr/fastqc.nf' params(run: true, outdir: params.outdir)
include { multiqc } from '../modules/crispr/multiqc.nf' params(run: true, outdir: params.outdir,
                                                   runtag : params.runtag)

// include workflow common to all input modes:
//include { run_from_irods_tsv } from './run_from_irods_tsv.nf'
//include { crispr } from './crispr.nf'

workflow {

    if (params.run_mode == "study_id") {
	imeta_study(Channel.from(params.study_id_mode.input_studies))
	samples_irods_tsv = imeta_study.out.irods_samples_tsv
	work_dir_to_remove = imeta_study.out.work_dir_to_remove }
    
    else if (params.run_mode == "id_run") {
             imeta_run(Channel.from(params.run_id_mode.input_id_runs))
   	     samples_irods_tsv = imeta_run.out.irods_samples_tsv
	     work_dir_to_remove = imeta_run.out.work_dir_to_remove }

    else if (params.run_mode == "csv_samples_id") {
	i1 = Channel.fromPath(params.csv_samples_id_mode.input_samples_csv)
	i2 = Channel.from(params.csv_samples_id_mode.input_samples_csv_column)
	imeta_samples_csv(i1,i2)
	samples_irods_tsv = imeta_samples_csv.out.irods_samples_tsv
	work_dir_to_remove = imeta_samples_csv.out.work_dir_to_remove }
    
    else if (params.run_mode == "google_spreadsheet") {
	i1 = Channel.from(params.google_spreadsheet_mode.input_gsheet_name)
	i2 = Channel.fromPath(params.google_spreadsheet_mode.input_google_creds)
	i3 = Channel.from(params.google_spreadsheet_mode.output_csv_name)
	gsheet_to_csv(i1,i2,i3)
	i4 = Channel.from(params.google_spreadsheet_mode.input_gsheet_column)
	imeta_samples_csv(gsheet_to_csv.out.samples_csv, i4)
	samples_irods_tsv = imeta_samples_csv.out.irods_samples_tsv
	work_dir_to_remove = imeta_samples_csv.out.work_dir_to_remove.mix(gsheet_to_csv.out.work_dir_to_remove) }

    // common to all input modes:

    
    //run_from_irods_tsv(samples_irods_tsv)
    //crispr()

    iget_study_cram(
        samples_irods_tsv
            .map{study_id, samples_tsv -> samples_tsv}
            .splitCsv(header: true, sep: '\t')
            .map{row->tuple(row.study_id, row.sample, row.object)}
            .filter { it[2] =~ /.cram$/ } // Need to check for bam too?
            .dump()
            .unique())

    crams_to_fastq(iget_study_cram.out.study_sample_cram.groupTuple(by: [0,1]))

    // store the number of reads in merged cram in output tables
    // lostcause has samples that did not pass the crams_to_fastq_min_reads input param, which is the minimum number of reads in merged cram file to try and convert to fastq.gz

    crams_to_fastq.out.lostcause
        .collectFile(name: "crams_to_fastq_lowreads.tsv",
                     newLine: false, sort: true, keepHeader: true,
                     storeDir:params.outdir)
    // numreads has all samples that pass min number of reads number of reads in merged cram file
    crams_to_fastq.out.numreads
        .collectFile(name: "crams_to_fastq_numreads.tsv",
                     newLine: false, sort: true, keepHeader: true,
                     storeDir:params.outdir)
    crams_to_fastq.out.fastq
                .collectFile(name: "samplename_to_fastq.csv",
                     newLine: false, sort: true, keepHeader: true,
                     storeDir:params.outdir)

   fastx_trimmer(crams_to_fastq.out.fastq_for_trim)

    fastx_trimmer.out
        .groupTuple(by: 0, sort: true)
        .map{ samplename, batchs, fastqs -> tuple( groupKey(samplename, batchs.size()), batchs, fastqs ) }
        .set{ch_samplename_fastqs_to_merge}

    merge_fastq_batches(ch_samplename_fastqs_to_merge)

    fastqc(fastx_trimmer.out
            .map{ samplename, batch, fastq -> tuple( samplename, fastq ) }
            .mix(merge_fastq_batches.out[0]))

    multiqc(fastqc.out.collect())

    merge_fastq_batches.out[0].view()

    merge_fastq_batches.out[0]
        .map{sample,fastq ->tuple("$sample",fastq)}
        .combine(ch_samplename_library
                 .map{sample,lib,keepg ->tuple("$sample",lib,keepg)}
                 , by: 0)
        .set{ch_samplename_fastq_library_includeG}

    count_crispr_reads(ch_samplename_fastq_library_includeG, ch_library_files.collect())

    collate_crispr_counts(
        count_crispr_reads.out[0]
            .map{ lib_csv,counts -> [ lib_csv.replaceAll(~/.csv/, ""), counts ] }
            .transpose()
            .groupTuple(by: 0, sort: true)
            .mix(count_crispr_reads.out[0].
                 map{lib,counts -> counts}.collect().map{a -> tuple("all_libs", a)})
    )

    count_crispr_reads.out[1].collectFile(name: 'mapping_percent.txt', newLine: true,
                                          storeDir: "${params.outdir}/", sort: true)



    // list work dirs to remove (because they are Irods searches, so need to always rerun on each NF run):
    // these are removed on workflow.onComplete if (params.on_complete_uncache_irods_search), see below.

    if (params.run_mode == "google_spreadsheet") {
	// combine all samples tables (google spreadsheet, irods + cellranger metadata, cellranger /lustre file paths),
	//   by joining on common sample column:
	// the resulting combined tables can be fed directly as input to the Vireo deconvolution pipeline or QC pipeline.
	join_gsheet_metadata(gsheet_to_csv.out.samples_csv,
			     run_from_irods_tsv.out.ch_cellranger_metadata_tsv,
			     run_from_irods_tsv.out.ch_file_paths_10x_tsv)
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
	def proc = "bash ${projectDir}/../bin/del_work_dirs_failed.sh ${workDir}".execute()
	def b = new StringBuffer()
	proc.consumeProcessErrorStream(b)
	log.info proc.text
	log.info b.toString() }
}

