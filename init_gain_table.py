import numpy as np

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
