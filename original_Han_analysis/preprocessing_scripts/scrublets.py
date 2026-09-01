#!/usr/bin/env python

import argparse as ap
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

# Set up the argument parser
argparser = ap.ArgumentParser(description="Run Scrublet doublet detection on one sample/lane")

argparser.add_argument( "-s", "--sample", dest = "sample_name", type = str,
                        help = "Sample/lane folder name, e.g. P0WT_lane1.",
                        default = "P0WT_lane1")
# Project root holds the per-sample folders. Defaults to the parent of this
# script's directory (i.e. scripts/.. == project root), so it is portable.
argparser.add_argument( "-p", "--project-dir", dest = "project_dir", type = str,
                        default = os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        help = "Project root containing the sample folders.")

args = argparser.parse_args()
sample_name = args.sample_name
project_dir = args.project_dir

input_dir = os.path.join(project_dir, sample_name, 'filtered_matrix', 'sensitivity_5')
output_dir = os.path.join(project_dir, sample_name)


counts_matrix = scipy.io.mmread(os.path.join(input_dir, 'matrix.mtx')).T.tocsc()
genes = np.array(scr.load_genes(os.path.join(input_dir, 'features.tsv'), delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))



scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.08)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

df = pd.DataFrame({
    'doublet_scores': doublet_scores,
    'predicted_doublets': predicted_doublets
})  

df.to_csv(os.path.join(output_dir, 'predicted_doublets.csv'), index=False)

# scrub.plot_histogram();
fig, ax = scrub.plot_histogram()
fig.savefig(os.path.join(output_dir, 'doublet_score_histogram.png'), dpi=300, bbox_inches='tight')

print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

# # Uncomment to run tSNE - slow
# print('Running tSNE...')
# scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

# # Uncomment to run force layout - slow
# print('Running ForceAtlas2...')
# scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5. n_iter=1000))
    
print('Done.')       

# fig, ax = scrub.plot_embedding('UMAP', order_points=True);
# fig.savefig(os.path.join(output_dir, 'Predicted_doublets_umap.png'), dpi=300, bbox_inches='tight')
                                                          