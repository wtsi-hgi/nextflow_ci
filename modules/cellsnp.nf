process cellsnp {
    tag "${samplename}"
    publishDir "${params.outdir}/cellsnp/", mode: 'copy', pattern: "cellsnp_${samplename}", overwrite: true

    input: 
    tuple val(samplename), path(bam_file), path(bai_file), path(barcodes_tsv_gz)
    file(region_vcf)
    
    output:
    tuple val(samplename), file("cellsnp_${samplename}"), emit: cellsnp_output_dir
    tuple val(samplename), file("cellsnp_${samplename}/cellSNP.cells.vcf.gz"), emit: samplename_cellsvcfgz

    script:
    """
umask 2 # make files group_writable

zcat ${barcodes_tsv_gz} > barcodes.txt
cellSNP -s ${bam_file} -b barcodes.txt -O cellsnp_${samplename} -R ${region_vcf} -p 20 --minMAF 0.1 --minCOUNT 60
    """
}
// https://github.com/single-cell-genetics/cellSNP
// Mode 1: pileup a list of SNPs for single cells in a big BAM/SAM file
// Require: a single BAM/SAM file, e.g., from cellranger, a VCF file for a list of common SNPs.
// This mode is recommended comparing to mode 2, if a list of common SNP is known, e.g., human (see Candidate SNPs)
// As shown in the above command line, we recommend filtering SNPs with <20UMIs or <10% minor alleles for downstream donor deconvolution, by adding --minMAF 0.1 --minCOUNT 20
