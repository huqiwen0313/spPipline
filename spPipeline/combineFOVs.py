import os, re, numpy as np, pandas as pd
from scipy.spatial import cKDTree

fov_pat = r"FOV(\d+)"

def removeOverlapRolonies(rolonyDf, x_col = 'x', y_col = 'y', removeRadius = 5.5):
    """ For each position, find those rolonies that are very close to other rolonies 
        in other positions and remove them.
        x_col and y_col are the names of the columns for x and y coordinates.
        removeRadius is in any unit that x_col and y_col are.
    """
    geneList = rolonyDf.target.unique()
    reducedRolonies = []
    for gene in geneList:
        thisGene_rolonies = rolonyDf.loc[rolonyDf.target == gene]
        for pos in sorted(rolonyDf['fov'].unique()):
            thisPos = thisGene_rolonies.loc[thisGene_rolonies['fov'] == pos]
            otherPos = thisGene_rolonies.loc[thisGene_rolonies['fov'] != pos]
            if (len(thisPos) <= 0 ) or (len(otherPos) <= 0 ):
                continue
            nnFinder = cKDTree(thisPos[[x_col, y_col]])
            nearestDists, nearestInds = nnFinder.query(otherPos[[x_col, y_col]], distance_upper_bound = removeRadius)
            toRemoveFromThisPos_index = thisPos.index[nearestInds[nearestDists < np.inf]]
            thisGene_rolonies = thisGene_rolonies.drop(toRemoveFromThisPos_index)
        reducedRolonies.append(thisGene_rolonies)
    return pd.concat(reducedRolonies) 


def filterByEmptyFraction(spot_df, cutoff):
    spot_df = spot_df.sort_values('distance')
    spot_df['isEmpty'] = spot_df['target'].str.startswith('Empty')
    spot_df['cum_empty'] = spot_df['isEmpty'].cumsum()
    spot_df['cum_empty_rate'] = spot_df['cum_empty'] / np.arange(1, spot_df.shape[0] + 1)
    spot_df_trimmed = spot_df.loc[spot_df['cum_empty_rate'] <= cutoff]
    return spot_df_trimmed, spot_df


def makeSpotTable(files_paths, emptyFractionCutoff, voxel_info):
    # Concatenating spots from all FOVs and converting the physical coordinates to pixels 
    allspots = []
    for file in files_paths: 
        thisSpots = pd.read_csv(file, index_col = 0)
        thisSpots['x'] = (round(thisSpots['xc'] / voxel_info['X'])).astype(int)
        thisSpots['y'] = (round(thisSpots['yc'] / voxel_info['Y'])).astype(int)
        thisSpots['z'] = (round(thisSpots['zc'] / voxel_info['Z'])).astype(int)
        thisSpots['fov'] = re.search(fov_pat, file).group()
        allspots.append(thisSpots)

    allspots = pd.concat(allspots, ignore_index=True)
    
    allspots['gene'] = allspots['target'].str.extract(r"^(.+)_")
    
    allspots = allspots.sort_values('distance')

    # Removing duplicate rolonies caused the overlapping regions of FOVs
    allspots_reduced = removeOverlapRolonies(allspots, x_col='x', y_col = 'y', removeRadius=5.5)

    # Keeping only spots with small distance to barcode so that `emptyFractionThresh` of spots are empty.
    allspots_trimmed, allspots_reduced = filterByEmptyFraction(allspots_reduced, cutoff=emptyFractionCutoff)

    return allspots_trimmed, allspots_reduced


def combine_fovs(decoding_dir, voxel, emptyFractionThresh=0.12):

    bcmags = [file for file in os.listdir(decoding_dir)
              if os.path.isdir(os.path.join(decoding_dir, file))
              and 'bcmag' in file]  # all the bcmags that were used for the experiment

    for bcmag in bcmags:
        print("filtering barcode magnitude: {}".format(bcmag))
        all_files = [os.path.join(decoding_dir, bcmag, file)
                 for file in os.listdir(os.path.join(decoding_dir, bcmag))
                 if re.search(fov_pat, file)]

        all_files.sort(key=lambda x: int(re.search(fov_pat, x).group(1)))
        filtered_spots, _ = makeSpotTable(all_files, emptyFractionThresh, voxel)
        filtered_spots.reset_index(drop=True).to_csv(os.path.join(decoding_dir, bcmag, 'all_spots_filtered.tsv'),
                                                   sep='\t')
