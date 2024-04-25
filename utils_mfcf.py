import numpy as np
import itertools as itertools

def id_from_set(set_matrix, an_element):
    partial = np.sum(np.isin(set_matrix, an_element), axis=1)
    idx = np.where(partial == np.sum(~np.isnan(an_element)))[0]
    return idx


def apply_threshold_v(ranked_values, ranked_seps, mincsize, threshold):

    ranked_values_thr = []
    ranked_seps_thr = []
    for i in range(len(ranked_values)):
        #rank_values_thr_row = []
        #rank_values_seps_row = []
        idx = np.where(ranked_values[i]<threshold)
        idx = (idx[0]>mincsize-2)*idx[0]
        ranked_values[i][idx]=0
        ranked_seps[i][idx]=np.nan

        ranked_values_thr.append(ranked_values[i])
        ranked_seps_thr.append(ranked_seps[i])  

    return np.array(ranked_values_thr),np.array(ranked_seps_thr)


def greedy_sortsep_v(vertices, sets, W):
    
    ranked_values = []
    ranked_seps = []
    #for i in range(len(the_vs)):
    for i in range(len(vertices)):
        val,sep =  greedy_sortsep(vertices[i], sets[i], W)
        ranked_values.append(val)
        ranked_seps.append(sep)

    return np.array(ranked_values), np.array(ranked_seps)


def greedy_sortsep(vtx, sep, W):
    sep_ranked = sep[~np.isnan(sep)]
    pad = len(sep) - len(sep_ranked)
    values, ranks = np.sort(W[vtx, sep_ranked])[::-1], np.argsort(W[vtx, sep_ranked])[::-1]
    sep_ranked = np.concatenate((sep_ranked[ranks], np.full(pad, np.nan)))
    values = np.concatenate((values, np.zeros(pad)))
    return values, sep_ranked


def j_LoGo(cov, cliques, separators):
    J = np.zeros((cov.shape[0],cov.shape[0]))
    C = cov#.to_numpy()
    for c in cliques:
        if np.isnan(c).any()==False:
            c = np.array(c,dtype=int)
            J[np.ix_(c, c)] += np.linalg.inv(C[np.ix_(c, c)])
    for s in separators:
        if np.isnan(s).any()==False:
            s = np.array(s,dtype=int)
            J[np.ix_(s, s)] -= np.linalg.inv(C[np.ix_(s, s)])

    return J

def init_gain_table(*args):
    # Optional parameters
    if len(args) == 2:
        MAX_CLIQUE_SIZE = args[0]
        INITIAL_GAIN_TABLE_SIZE = args[1]
    elif len(args) == 1:
        MAX_CLIQUE_SIZE = args[0]
        INITIAL_GAIN_TABLE_SIZE = 100
    else:
        MAX_CLIQUE_SIZE = 6
        INITIAL_GAIN_TABLE_SIZE = 100

    MAX_SEPARATOR_SIZE = MAX_CLIQUE_SIZE - 1

    GT = {
        'MAX_CLIQUE_SIZE': MAX_CLIQUE_SIZE,
        'INITIAL_GAIN_TABLE_SIZE': INITIAL_GAIN_TABLE_SIZE,
        'MAX_SEPARATOR_SIZE': MAX_SEPARATOR_SIZE- 1,
        'rowid': np.full(INITIAL_GAIN_TABLE_SIZE, np.nan),
        'cliques': np.full((INITIAL_GAIN_TABLE_SIZE, MAX_CLIQUE_SIZE), np.nan),
        'separators': np.full((INITIAL_GAIN_TABLE_SIZE, MAX_SEPARATOR_SIZE), np.nan),
        'coordination_num': np.full(INITIAL_GAIN_TABLE_SIZE, np.nan),
        'gains': np.full(INITIAL_GAIN_TABLE_SIZE, np.nan),
        'nodes': np.full(INITIAL_GAIN_TABLE_SIZE, np.nan),
        'tot_records': 0
    }

    return GT

