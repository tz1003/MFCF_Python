import numpy as np
from itertools import combinations
import itertools as itertools
import pandas as pd
from scipy.sparse import lil_matrix
from gain_function import gf_sumsquares_gen
from utils_mfcf import *

def MFCF_Forest(X,ct_control,gain_function):

    # analyse X
    X = np.array(X)
    sizeX = X.shape
    isSquare = sizeX[0] == sizeX[1]
    p = sizeX[1]
    max_cliques = p - 1

    # Preallocate output
    cliques = np.full((max_cliques, ct_control['max_clique_size']), np.nan)
    separators = np.full((max_cliques - 1, ct_control['max_clique_size'] - 1), np.nan)
    peo = []
    tree = lil_matrix((max_cliques, max_cliques))
    GT = init_gain_table(ct_control['max_clique_size'])

    outstanding_nodes = np.arange(p)
    clique_no = 0
    sep_no = 0

    # If matrix is not square, make a kernel matrix
    if not isSquare:
        kernel_func = lambda r, c: gain_function(X, r, c, ct_control) #########working to here
        ij = np.array(list(combinations(range(1, p + 1), 2)))
        vals = np.array([kernel_func(ij[:,0], ij[:,1])])
        K = np.zeros((p, p))
        idx = np.ravel_multi_index((ij[:, 0] - 1, ij[:, 1] - 1), K.shape)
        K.flat[idx] = vals
        K = K + K.T
        sums = np.sum(K * (K > np.mean(K)), axis=1)
    else:
        sums = np.sum(X * (X > np.mean(X)), axis=1)


    # Get first ct_control['min_clique_size'] nodes as the first clique
    j = np.argsort(sums)[::-1]
    first_clique = j[:ct_control['min_clique_size']]
    clique_no = 1
    cliques[clique_no - 1, :len(first_clique)] = first_clique
    # Remove first clique from outstanding_nodes
    outstanding_nodes = outstanding_nodes[~np.isin(outstanding_nodes, first_clique)]
    peo = first_clique.copy()
    # Calculate gains for this clique
    nodes, gains, seps = gain_function(X, first_clique, outstanding_nodes, ct_control)
    new_gains = len(gains)
    from_idx, to_idx = GT['tot_records'], GT['tot_records'] + new_gains
    GT['tot_records'] = to_idx
    idx = np.arange(from_idx, to_idx)
    clq_len = len(first_clique)


    if np.isnan(GT['cliques'][:,0]).sum()<len(idx):
        diff = (len(idx)-np.isnan(GT['cliques'][:,0]).sum())*2
        adding_item_clique = np.full((diff, ct_control['max_clique_size']), np.nan)
        GT['cliques']=np.concatenate((GT['cliques'], adding_item_clique), axis=0)
        adding_item_gains = np.array([np.nan]*diff)
        GT['gains']=np.concatenate((GT['gains'], adding_item_gains), axis=0)
        adding_item_separators = np.full((diff, ct_control['max_clique_size']-1), np.nan)
        GT['separators']=np.concatenate((GT['separators'], adding_item_separators), axis=0)


    GT['cliques'][idx, :clq_len] = np.tile(first_clique, (new_gains, 1))
    if clq_len < ct_control['max_clique_size']:
        GT['cliques'][idx, clq_len:] = np.nan
    GT['gains'][idx] = gains # this result doesnt seem same
    GT['separators'][idx, :] = seps
    GT['nodes'] = nodes


    while len(outstanding_nodes) > 0:

        the_gain = np.nanmax(GT['gains'])

        # Case: no gain, add an isolated clique
        if np.isnan(the_gain):
            the_node = outstanding_nodes[0]
            the_sep = np.array([], dtype=int)
            parent_clique = np.array([], dtype=int)
            parent_clique_id = np.nan
        else:
            idx = np.nanargmax(GT['gains'])
            the_node = int(GT['nodes'][idx])
            the_sep = GT['separators'][idx, :]
            the_sep = np.array(the_sep[~np.isnan(the_sep)],dtype=int)
            
            parent_clique = GT['cliques'][idx, :]
            parent_clique_id = id_from_set(cliques, parent_clique)

        new_clique = np.concatenate((the_sep, [the_node]),dtype=int)
        outstanding_nodes = outstanding_nodes[outstanding_nodes != the_node]
        peo = np.concatenate((peo, [the_node]))

        # Track if it is a clique extension or a new clique
        clique_extension = 0

        # In case of no gain, add an isolated clique
        if np.isnan(the_gain):
            clique_no += 1
            clique_to_update = clique_no
            cliques[clique_to_update - 1, :len(new_clique)] = new_clique
        # Case: new clique with an existing intersection
        elif len(new_clique) <= np.sum(~np.isnan(parent_clique)):
            clique_no += 1
            clique_to_update = clique_no
            sep_no += 1
            cliques[clique_to_update - 1, :len(new_clique)] = new_clique
            separators[sep_no - 1, :len(the_sep)] = the_sep
            tree[clique_to_update - 1, parent_clique_id] = 1
        # Case: extension of an existing clique
        else:
            clique_to_update = parent_clique_id
            cliques[clique_to_update, :len(new_clique)] = new_clique
            clique_extension = 1
            # Parent clique is the same as the clique without the last extension
            old_clique_idx = id_from_set(GT['cliques'], parent_clique)

        # Don't update the gain table when nodes are finished
        if len(outstanding_nodes) == 0:
            break

        # Finally, update the gain table
        nodes, gains, seps = gain_function(X, new_clique, outstanding_nodes, ct_control)
        new_gains = len(gains)
        from_idx, to_idx = GT['tot_records'], GT['tot_records'] + new_gains
        GT['tot_records'] = to_idx

        idx = np.arange(from_idx, to_idx)
        new_clique = new_clique[~np.isnan(new_clique)]
        clq_len = len(new_clique)

        if np.isnan(GT['cliques'][:,0]).sum()<len(idx):
            diff = len(idx)-np.isnan(GT['cliques'][:,0]).sum()
            adding_item_clique = np.full((diff, ct_control['max_clique_size']), np.nan)
            GT['cliques']=np.concatenate((GT['cliques'], adding_item_clique), axis=0)
            adding_item_gains = np.array([np.nan]*diff)
            GT['gains']=np.concatenate((GT['gains'], adding_item_gains), axis=0)
            adding_item_separators = np.full((diff, ct_control['max_clique_size']-1), np.nan)
            GT['separators']=np.concatenate((GT['separators'], adding_item_separators), axis=0)


        GT['cliques'][idx, :clq_len] = np.tile(new_clique, (new_gains, 1))
        if clq_len < ct_control['max_clique_size']:
            GT['cliques'][idx, clq_len:] = np.nan
        GT['gains'][idx] = gains
        GT['separators'][idx, :] = seps
        GT['nodes'] = np.concatenate((GT['nodes'], nodes), axis=0)

        # Remove from the gain table where the node is the new one
        idx = np.where(GT['nodes'] == the_node)[0]
        GT['gains'][idx] = np.nan

        # If the clique was expanded, remove the records with the old clique from the gain table
        if clique_extension == 1:
            GT['gains'][old_clique_idx] = np.nan

        if len(the_sep)!=0:
            idx = id_from_set(separators, the_sep)
            if len(idx)>=ct_control['coordination_num']:
                idx = id_from_set(GT['separators'], the_sep)
                GT['gains'][idx] = np.nan

        # If drop sep, remove also the separator just used
        if ct_control['drop_sep'] == True:
            idx = id_from_set(GT['separators'], the_sep)
            GT['gains'][idx] = np.nan

    peo = np.flipud(peo)

    return cliques, separators, peo, tree
