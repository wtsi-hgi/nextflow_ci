#!/usr/bin/env python3
__author__ = 'Guillaume Noell'
__date__ = '2021-01-15'
__version__ = '0.0.1'
# for help, run: python3 plot_donor_ncells.py --help

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
@click.option('--sample_donor_summary_tsv', required=True, type=click.Path(exists=True),
              help='path to multi-samples donor_ids.tsv file, which is a concatenation of multiple donor_ids.tsv with an additional column for sample (10x experiment ID)')

@click.option('-d','--plotnine_dpi', default=100,
              type=click.IntRange(1, 1000, clamp=True), show_default=True,
              help='DPI pdf plots resolution for plotnine plots. Integer in range 1 to 1000')


def plot_donor_ncells(sample_donor_summary_tsv, plotnine_dpi):
    """plot_donor_ncells main script"""
    logging.info('running plot_donor_ncells() function..')

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
    plot_donor_ncells()
