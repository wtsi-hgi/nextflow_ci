process subset_genotype {
    tag "${samplename}"
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.subset_genotype.copy_mode}", pattern: "${samplename}.subset.vcf.gz"
    
    when: 
    params.vireo_with_genotype.subset_genotype.run

    when:
    params.run 

    input:
    set val(samplename), file(cellsnp_vcf), file(donor_vcf), file(sample_subset_file)
    
    output:
    tuple val(samplename), file("${samplename}.subset.vcf.gz"), emit: samplename_subsetvcf

    script:
    """
tabix -p vcf ${donor_vcf} 
# tabix -p vcf ${cellsnp_vcf}
bcftools view ${donor_vcf} -R ${cellsnp_vcf} -S ${sample_subset_file} -Oz -o ${samplename}.subset.vcf.gz
rm ${donor_vcf}.tbi
# rm ${cellsnp_vcf}.tbi
    """

}
