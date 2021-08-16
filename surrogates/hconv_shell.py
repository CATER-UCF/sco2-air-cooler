"""
A simple linear fit is sufficient:

hconv = 5.00519828 + 0.04982918 * temp
"""
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

npr = 20
nt = 20

df = pd.read_csv(f'./data/DOE_air_p{npr}_t{nt}.csv')

hconv = np.array(df['hconv'])
temp = np.array(df['temperature'])


def poly(x, a0, a1):
    return a0 + a1 * x #+ a2 * x ** 2


popt, _ = optimize.curve_fit(poly, temp, hconv)

y_fit = poly(temp, *popt)
print(popt)

fig, ax = plt.subplots(2)
ax[0].plot(temp, hconv, label='hconv')
ax[0].plot(temp, y_fit, '.', label='fit')
ax[0].legend()

ax[1].plot(temp, y_fit - hconv, '.', label='error')
ax[1].legend()
plt.show()
