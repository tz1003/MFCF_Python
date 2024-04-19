# MFCF_Python

repository for MFCF Python version

The package is still under development, more examples/comments will be added soon


## Progress

✔ Main function

✔ Examples

Function *drop_sep*

Bug on *coord_num*

Visualisation of the networks



## Parameters
**max_clique_size**: maximum clique size

**min_clique_size**:minimum clique size

**threshold**: only edges with score above the threshold will appear in the final result

**coordination_num**: number below which each nodes are allowed to connect with

**drop_sep**:use any separator only once (suggest no)


## Usage Example
```python
from MFCF_main import MFCF_Forest
from gain_function import gf_sumsquares_gen
import numpy as np
import pandas as pd

def generate_random_df(num_rows, num_columns):
  data = np.random.randint(0, 100, size=(num_rows, num_columns))
  df = pd.DataFrame(data, columns=['col_{}'.format(i) for i in range(num_columns)])
  return df

# fill in ct_control
ct_control = {
    'max_clique_size': 2,
    'min_clique_size': 1,
    'threshold': 0.00,
    'coordination_num':np.inf,
    'drop_sep': False
}

df = generate_random_df(100, 100)
cliques, separators, peo, tree = MFCF_Forest(df,ct_control,gf_sumsquares_gen)
```
More examples could be found in the folder *examples*

## Feedback

For any feedback, kindly create an issue
