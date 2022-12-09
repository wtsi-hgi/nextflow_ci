#!/bin/bash
#  b build_star_index_2.7.8a.sh 12 20 normal
star_version=2.7.8a
singularity_cache=/lustre/scratch125/humgen/resources/containers
# from https://github.com/wtsi-hgi/nextflow_rna_seq_container : 
image=rna_seq_1.0.img

n_threads=12
star_index_dir=star${star_version}_index_Homo_sapiens.GRCh38.99_75bp
mkdir ${star_index_dir}
singularity exec --bind /lustre --bind /tmp ${singularity_cache}/${image} /bin/bash -c \
	    "cd $PWD; STAR \
        --runMode genomeGenerate \
        --runThreadN ${n_threads} \
        --sjdbGTFfile Homo_sapiens.GRCh38.99.gtf \
        --sjdbOverhang 75 \
        --genomeDir ${star_index_dir} \
        --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa"

star_index_dir=star${star_version}_index_Homo_sapiens.GRCh38.99_100bp
mkdir ${star_index_dir}
singularity exec --bind /lustre --bind /tmp ${singularity_cache}/${image} /bin/bash -c \
	    "cd $PWD; STAR \
        --runMode genomeGenerate \
        --runThreadN ${n_threads} \
        --sjdbGTFfile Homo_sapiens.GRCh38.99.gtf \
        --sjdbOverhang 100 \
        --genomeDir ${star_index_dir} \
        --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa"
