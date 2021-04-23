nextflow.enable.dsl=2

// include RNAseq Nextflow modules/tasks:
include { iget_cram } from '../modules/irods.nf'
include { crams_to_fastq_gz } from '../modules/crams_to_fastq.nf'
include { baton_study_id } from '../modules/baton.nf'
include { mbv } from '../modules/mbv.nf'
include { get_egan_id } from '../modules/get_egan_id.nf'

workflow get_fastq {
    
    take:
    ch_input_study_ids // input Irods study IDs to get CRAM files from
    ch_input_fastq_csv // input fastq csv file
    
    main:
    log.info "\nrunning get_fastq() sub-workflow to get fastq files, depending on input mode (from Irods CRAM files or from list of fastq files)..."

    if (params.input_mode == "from_study_id") {
	log.info "params.input_mode == \"from_study_id\"" 
	// Get list of study cram files from Irods:
	baton_study_id(ch_input_study_ids)
	
	// Download cram files from Irods:
	to_iget = baton_study_id.out.samples_noduplicates_tsv
	    .map{a,b -> b}
	    .splitCsv(header: true, sep: '\t')
	    .map{row->tuple(row.sample, row.sample_supplier_name, row.study_id)}
	    .map{a,b,c-> tuple(a,c)}
	iget_cram(to_iget)
	
	// Optionally, extract egan IDs from cram files:
	get_egan_id(iget_cram.out[0])
	get_egan_id.out.samplename_egan_id_csv
	    .map { samplename,cram,csv_file -> csv_file }
	    .splitCsv(header: true, sep: ',')
	    .map { row -> "${row.samplename},${row.egan_id}"}
	    .collectFile(name: 'samplename_egan_id.csv', newLine: true,
			 seed: "samplename,egan_id",
			 storeDir: "${params.outdir}/", sort: true)
	    .set{ch_samplename_egan_id_csv}

	// Convert CRAM files to fastq:
	crams_to_fastq_gz(iget_cram.out[0])
	crams_to_fastq_gz.out[0]
	    .map{ samplename, fastq1, fastq2 -> tuple( samplename, tuple(fastq1, fastq2) ) }
            .set{ ch_samplename_crams }
	
    } else if (params.input_mode == "from_fastq_csv") {
	log.info "params.input_mode == \"from_fastq_csv\"" 
	Channel
	    .fromPath(params.from_fastq_mode.fastq_csv)
	    .splitCsv(header:false)
	    .map{ row-> tuple(row[0], tuple(file(row[1]), file(row[2]))) }
	    .set { ch_samplename_crams }
    } else if (params.run_mode != "from_fastq_csv" & params.run_mode == "from_study_id") {
        log.info "error: input param run_mode must be set to either \"from_study_id\" or \"from_fastq_csv\""
        exit 1
    }
    
    emit:
    ch_samplename_crams // channel of tuple(samplename,tuple(fastq1/fastq2)) for paired end reads.
    
}

// possible additional input mode to include,
//// from cram files:

    ////Channel.fromPath('/lustre/scratch115/projects/interval_rna/inputs/*.cram').
    ////map{ it -> [ it.toString().replaceAll(~/.*\/(.*).cram/, "\$1"), it ] }.
    ////groupTuple(sort: true). //take(4).
    ////set{ch_cram_files}
    ////crams_to_fastq_gz(ch_cram_files)
    ////
    //crams_to_fastq_gz.out[0]
    //	.map{ samplename, fastq1, fastq2 -> tuple( samplename, tuple(fastq1, fastq2) ) }
    //	.set{ch_samplename_crams}

//// from cram files:
////Channel.fromPath('/lustre/scratch115/projects/interval_rna/inputs/*.cram').
////map{ it -> [ it.toString().replaceAll(~/.*\/(.*).cram/, "\$1"), it ] }.
////groupTuple(sort: true). //take(4).
////set{ch_cram_files}
////crams_to_fastq_gz(ch_cram_files)
