import numpy as np
import itertools as itertools
import pandas as pd
from utils_mfcf import *

def gf_sumsquares_gen(M, clq, v, ct_control):
    """
    M similarity matrix
    v vector of outstanding nodes
    clq clique
    ct_control parameters for the clique expansion algorithm
    """
    nodes = []
    gains = []
    seps = []

    clq = tuple(clq)
    csz = len(clq)
    vn = len(v)
    W = M*M

    if 'threshold' in ct_control:
        threshold = ct_control['threshold']
    else:
        threshold = 0.0
        
    if 'min_clique_size' in ct_control:
        min_clique_size = ct_control['min_clique_size']
    else:
        min_clique_size = 4

    if 'max_clique_size' in ct_control:
        max_clique_size = ct_control['max_clique_size']
    else:
        max_clique_size = 4

    if 'cachesize' in ct_control:
        cachesize = ct_control['cachesize']
    else:
        cachesize = min(4, min_clique_size)
        
    if csz < max_clique_size:
        facets = list()
        facets.append(clq)
    else:
        facets = list(itertools.combinations(clq, csz-1))

    block_rows = len(facets)
    ncols = len(facets[0])

    the_vs = np.sort(np.tile(v, block_rows)) # nodes in order
    the_fs = np.array(facets * vn,dtype=int) # facets as they are generated
    #repeated_matrix = np.tile(facets, (vn,1))


    ranked_values, ranked_seps = greedy_sortsep_v(the_vs, the_fs, W)
    ranked_values_thr, ranked_seps_thr =  apply_threshold_v(ranked_values, ranked_seps, min_clique_size, threshold)

    for i in range(len(the_vs)):
        gains.append(ranked_values_thr[i].sum())

    gains = np.array(gains)

    selector = np.tile(np.arange(1, block_rows + 1), vn)
    # Assuming the_vs, gains, ranked_seps_thr are numpy arrays
    the_table = np.column_stack((the_vs, gains.T, ranked_seps_thr))
    # Convert numpy array to pandas DataFrame for sorting
    df = pd.DataFrame(the_table)
    # Sort by first column ascending and second column descending
    df.sort_values(by=[0, 1], ascending=[True, False], inplace=True)
    # Assuming selector and cachesize are defined
    df = df[selector <= cachesize]
    #the_table = np.array(df)
    #nodes = the_table[:, 0]
    #gains = the_table[:, 1]
    #seps = np.column_stack([the_table[:, 2:], np.full((len(nodes), max_clique_size - 1 - ncols), np.nan)])

    # Assuming the_table is a pandas DataFrame
    the_table = df
    nodes = the_table.iloc[:, 0]
    gains = the_table.iloc[:, 1]
    # Assuming maxcsize and ncols are defined
    seps = pd.concat([the_table.iloc[:, 2:].reset_index(drop=True), pd.DataFrame(np.nan, index=np.arange(len(nodes)), columns=np.arange(max_clique_size - 1 - ncols))], axis=1)

    return np.array(nodes,dtype=int), gains, np.array(seps)#np.array(seps,dtype=int)