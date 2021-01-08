import os
from os import path
from spPipeline.code_lib.Segmentation_201019 import Segmentor2D
from spPipeline.code_lib.Assignment_201020 import *
from skimage.io import imread
import numpy as np, pandas as pd


def mask2centroid(maskImg):
    centroids = []
    for i in range(1, maskImg.max() + 1):
        xs, ys = np.where(maskImg == i)
        xc, yc = xs.mean().astype(int), ys.mean().astype(int)
        centroids.append((xc, yc))
    return np.array(centroids)


def segmentation():
    nuc_path = '../2_Registered/stitched/MIP_7_DRAQ5_ch00.tif'  # parameters here are subject to change

saving_path = '../4_CellAssignment'

bcmag = 'bcmag2.0'  # parameters here are subject to change

if not path.exists(saving_path):
    os.makedirs(saving_path)

spot_file = '../3_Decoded/output_Starfish/{}/all_spots_filtered.tsv'.format(bcmag)

nuc_img = imread(nuc_path)

# segmenting the nuclear image
segmentor = Segmentor2D()
mask = segmentor.segment([nuc_img], diameters=40,
                         out_files=[path.join(saving_path, 'segmentation_mask.npy')])[0]

# Rolony assignment
spot_df = pd.read_csv(spot_file, index_col=0, sep='\t')
assigner = RolonyAssigner(nucleiImg=mask, rolonyDf=spot_df, axes=['y', 'x'])
labels, dists = assigner.getResults()

spot_df['nucleus_label'] = labels
spot_df['dist2nucleus'] = np.round(dists, 2)
spot_df = spot_df.sort_values('nucleus_label', ignore_index=True)
spot_df.to_csv(path.join(saving_path, 'spots_assigned.tsv'), sep='\t', index=False, float_format='%.3f')

# plotting assigned rolonies
fig = plt.figure(figsize=(int(mask.shape[0] / 200), int(mask.shape[1] / 200)))
ax = fig.gca()
plotRolonies2d(spot_df, mask, coords=['x', 'y'], ax=ax, backgroudImg=nuc_img)
fig.savefig(path.join(saving_path, 'assigned_rolonies.png'),
            transparent=True, dpi=400, bbox_inches='tight')

# finding the nuclei centroids
centroids = mask2centroid(mask)
centroid_df = pd.DataFrame({'nucleus_label': np.arange(1, mask.max() + 1),
                            'centroid_x': centroids[:, 0], 'centroid_y': centroids[:, 1]})
centroid_df.to_csv(path.join(saving_path, 'nuclei_locations.tsv'), sep='\t', index=False)

# plotting the nuclei with their label
fig = plt.figure(figsize=(int(mask.shape[0] / 200), int(mask.shape[1] / 200)))
ax = fig.gca()
ax.imshow(nuc_img, cmap='gray')
ax.scatter(centroids[:, 1], centroids[:, 0], s=1, c='red')
for i in range(centroids.shape[0]):
    ax.text(centroids[i, 1], centroids[i, 0], str(i), fontsize=5, c='orange')
fig.savefig(path.join(saving_path, 'nuclei_map.png'),
            transparent=True, dpi=400, bbox_inches='tight')

# Making the cell by gene matrix
nuc_gene_df = spot_df[['nucleus_label', 'gene']].groupby(by=['nucleus_label', 'gene'], as_index=False).size()
nuc_gene_df = nuc_gene_df.reset_index().pivot(index='nucleus_label', columns='gene', values='size').fillna(0).astype(
    int)
nuc_gene_df.to_csv(path.join(saving_path, 'nucleus-by-gene.tsv'), sep='\t')
