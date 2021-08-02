"""
Surrogate model using a multivariate, piecewise linear function from Pyomo.
This is not yielding a good model thus far. Also, I'm not sure how to actually
implement this as a constraint in a unit model.

Submitted an issue here:
https://stackoverflow.com/questions/68626937/what-is-the-correct-usage-of-transformedpiecewiselinearfunctionnd-in-pyomo
"""
from pyomo.core.kernel.piecewise_library.transforms_nd import (
    PiecewiseLinearFunctionND,
    TransformedPiecewiseLinearFunctionND
)
from pyomo.environ import Var, Constraint
from scipy.spatial import Delaunay
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


npr = 30
nt = 99

df_star = pd.read_csv(f'./data/p{npr}_t_star.csv')
p_space = np.array(df_star['pressure']).reshape(-1, 1)
t_star = np.array(df_star['T*'])

df = pd.read_csv(f'./data/DOE_p{npr}_t{nt}.csv')

df = df[df['temperature'] > 300]
X = np.array(df[['temperature', 'pressure']])
Y = np.array(df['hconv'])
hconv_surr = PiecewiseLinearFunctionND(Delaunay(X), Y)

fig, ax = plt.subplots(3)

for i in range(npr):
    pi = p_space[i]
    df_tmp = df[abs(df['pressure'] - pi) < 0.01].sort_values('temperature')
    ti = np.array(df_tmp['temperature'])
    hi = np.array(df_tmp['hconv'])

    h_surr = [hconv_surr(np.array([t, pi], dtype=object)) for t in ti]

    ax[0].plot(ti, hi, label=f'p={pi}')
    ax[1].plot(ti, h_surr, label=f'p={pi}')

plt.show()
