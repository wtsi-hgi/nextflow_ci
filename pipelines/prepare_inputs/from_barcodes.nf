nextflow.enable.dsl=2

workflow from_barcodes {
    take: channel_input_data_table
    main:
    log.info "running workflow from_barcodes() ..."
    
    channel_input_data_table
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
	.map{row->tuple(row.experiment_id, row.n_pooled)}
	.set{ch_experiment_npooled}
    
    channel_input_data_table
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
	.map{row->tuple(row.experiment_id, row.data_path_bam_file ,row.data_path_barcodes)}
	.set{pre_ch_experiment_bam_barcodes}

    channel_input_data_table
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
	.map{row->tuple(row.experiment_id, row.data_path_filt_h5)}
	.set{pre_ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'

    // if params.replace_in_path set to true:
    // used if path to input files are mounted differently on the file-system
    // (e.g. if /lustre is mounted in Openstack on a different absolute path in nextflow  worker nodes)
    if (params.replace_in_path) {
	pre_ch_experiment_filth5
	    .map{experiment, path -> tuple(experiment, path.replaceFirst(/${params.replace_in_path_from}/,
									 params.replace_in_path_to))}
	    .set{ch_experiment_filth5}

	pre_ch_experiment_bam_barcodes
	    .map{experiment, pathbam, pathbarcodes -> tuple(experiment, 
							    pathbam.replaceFirst(/${params.replace_in_path_from}/,
										 params.replace_in_path_to),
							    pathbam.replaceFirst(/${params.replace_in_path_from}/,
										 params.replace_in_path_to).replaceFirst(/$/, ".bai"),
							    pathbarcodes.replaceFirst(/${params.replace_in_path_from}/,
										      params.replace_in_path_to))}
	    .set{ch_experiment_bam_bai_barcodes}
    } else {
	pre_ch_experiment_filth5
	    .set{ch_experiment_filth5}
	    
	ch_experiment_bam_barcodes
	    .map { a,b,c,d -> tuple(a, file(b), file("${b}.bai"), file(c))}
	    .set {ch_experiment_bam_bai_barcodes}
    }
    
    emit:
    ch_experiment_bam_bai_barcodes
    ch_experiment_npooled
    ch_experiment_filth5
}
