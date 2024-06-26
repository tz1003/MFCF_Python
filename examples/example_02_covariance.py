from MFCF_main import MFCF_Forest
from gain_function import gf_sumsquares_gen
from utils_mfcf import j_LoGo
import numpy as np
import pandas as pd

def generate_random_df(num_rows, num_columns):
  data = np.random.randint(0, 100, size=(num_rows, num_columns))
  df = pd.DataFrame(data, columns=['col_{}'.format(i) for i in range(num_columns)])
  return df

# set up ct_control
ct_control = {
    'max_clique_size': 2,
    'min_clique_size': 1,
    'threshold': 0.00,
    'coordination_num':np.inf,
    'drop_sep': False
}

# generate random df
df = generate_random_df(100, 100)
corr = np.square(df.corr())
cov = df.cov()
W = np.array(corr**2)
# apply mfcf
cliques, separators, peo, tree = MFCF_Forest(W,ct_control,gf_sumsquares_gen)
# apply J_LoGo
J = j_LoGo(np.array(cov), cliques, separators)
adj_matrix = pd.DataFrame(J, columns=['col_{}'.format(i) for i in range(len(J[0]))])   
