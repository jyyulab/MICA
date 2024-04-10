#!/usr/bin/env python3

import math

#%%
math.exp(-2)

#%%
math.exp(-1.5)

#%%
math.exp(-1)

#%%
math.exp(-0.5)

#%%
math.exp(0)

#%%
math.exp(0.5)

#%%
math.exp(1.0)

#%%
math.exp(2.0)

#%%
import numpy as np
reso_lst = []
for reso in list(np.arange(-2.0, 2.0 + 0.1, 0.2)):
    reso_lst.append(math.exp(reso))
print(reso_lst)
