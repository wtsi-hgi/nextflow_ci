#!/usr/bin/env python3
__author__ = 'Guillaume Noell'
__date__ = '2021-01-15'
__version__ = '0.0.1'
# for help, run: python3 split_h5ad_per_donor.py --help
# command with example inputs arguments: python3 split_h5ad_per_donor.py --samplename ukbb_scrna9479582 --vireo_donor_ids_tsv /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/nf_ci_scrna_deconv/results/vireo/vireo_ukbb_scrna9479582/donor_ids.tsv --filtered_matrix_h5 /lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/fetch_Submission_Data_Pilot_UKB/nextflow_ci/pipelines/../../results/iget_study_cellranger/5933/ukbb_scrna9479582/cellranger_ukbb_scrna9479582/filtered_feature_bc_matrix.h5

# import python libraries:
# on farm5, these libraries are installed in a conda environment: conda activate nextflow
import logging
import click
import sys 
import argparse
import os
import csv
import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import plotnine as plt9
from plotnine.ggplot import save_as_pdf_pages

# CLI arguments:
@click.command()

# required arguments:
@click.option('--vireo_donor_ids_tsv', required=True, type=click.Path(exists=True),
              help='path to donor_ids.tsv file, which is an output file of Vireo')

@click.option('--filtered_matrix_h5', required=True, type=click.Path(exists=True),
              help='path to filtered cells h5 data file created by cellranger count, which was deconvoluted by Vireo')

@click.option('--samplename', required=True, type=str,
              help='sample name of cellranger experiment deconvoluted by Vireo')

# arguments that are optional because they have default value:
@click.option('-o','--output_dir', default='./vireo_deconv_out', show_default=True, type=str,
              help='output directory for script output plots and files')

@click.option('-m','--print_modules_version', default=False, show_default=True, type=bool,
              help='True or False: whetherer to print version of all python module to file.')

@click.option('-p','--plot_n_cells_per_vireo_donor', default=True, show_default=True,  type=bool,
              help='True or False: whetherer to plot number of cells per deconvoluted donor to pdf in --output_dir')

@click.option('-w','--write_donor_level_filtered_cells_h5', default=True, show_default=True, type=bool,
              help='True or False: whetherer to write donor level scanpy hdf5 objects to dir --output_dir')

@click.option('-c','--anndata_compression_level', default=6,
              type=click.IntRange(1, 9, clamp=True), show_default=True,
              help='Gzip compression level for scanpy write of AnnData hdf5 objects. Integer in range 1 to 9')

@click.option('-d','--plotnine_dpi', default=100,
              type=click.IntRange(1, 1000, clamp=True), show_default=True,
              help='DPI pdf plots resolution for plotnine plots. Integer in range 1 to 1000')


def split_h5ad_per_donor(vireo_donor_ids_tsv, filtered_matrix_h5, samplename,
                         output_dir, print_modules_version, plot_n_cells_per_vireo_donor,
                         write_donor_level_filtered_cells_h5, plotnine_dpi,
                         anndata_compression_level):
    """split_h5ad_donor main script"""
    logging.info('running split_h5ad_per_donor() function..')

    # Set seed for reproducibility
    seed_value = 0
    # 0. Set `PYTHONHASHSEED` environment variable at a fixed value
    # os.environ['PYTHONHASHSEED']=str(seed_value)
    # 1. Set `python` built-in pseudo-random generator at a fixed value
    random.seed(seed_value)
    # 2. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)
    sns.set(style='whitegrid')
    # Set the default dpi
    plt9.options.dpi = plotnine_dpi   

if __name__ == '__main__':
    # set logging level and handler:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[logging.StreamHandler()]) # logging.FileHandler("debug.log"),
    split_h5ad_per_donor()

