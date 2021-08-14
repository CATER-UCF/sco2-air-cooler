"""
Surrogate model of dP / l using density:

dP / l =  2.19442286e+01 -2.68411214e+05 / rho +
         -7.17992324e-02 * rho + 3.93218678e-05 * rho ** 2

This produces an excellent fit with error <= 0.15%.
"""
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

npr = 30
nt = 99


df_star = pd.read_csv(f'./data/p{npr}_t_star.csv')
p_space = np.array(df_star['pressure']).reshape(-1, 1)

df = pd.read_csv(f'./data/DOE_p{npr}_t{nt}.csv')
df = df[df['temperature'] > 299]


def kernel(x, a0, a1, a2, a3):
    """
    dP / L = a0 + a1 / rho + a2 * rho + a3 * rho ** 2
    """
    return a0 + a1 / x[0] + a2 * x[0] + a3 * x[0] ** 2


def get_model(df):
    # Initial guesses
    a0 = -2000
    a1 = 1
    a2 = 0.1
    a3 = 0.01
    p_init = (a0, a1, a2, a3)

    XX = np.empty((2, len(df)))
    XX[0, :] = np.array(df['rho'])
    XX[1, :] = np.array(df['pressure'])
    temps =  np.array(df['temperature'])
    Y = np.array(df['dP_over_l'])
    popt, _ = optimize.curve_fit(kernel, XX, Y, p0=p_init)
    return popt, XX, Y, temps


popt_generic, _, _, _ = get_model(df)
print(popt_generic)

fig, ax = plt.subplots(2)

popts = []
errs = np.array([])
errs_generic = np.array([])

for p in p_space:
    df_tmp = df[abs(df['pressure'] - p[0]) < 0.01]
    popt, XX, Y, temps = get_model(df_tmp)
    Y_model = kernel(XX, *popt)
    Y_generic = kernel(XX, *popt_generic)
    err = (Y_model - Y) / Y
    err_generic = (Y_generic - Y) / Y
    ax[0].plot(temps, err, '.')
    ax[1].plot(temps, err_generic, '.')
    popts.append(popt)
    errs = np.append(errs, err)
    errs_generic = np.append(errs_generic, err_generic)


fig, ax = plt.subplots(2)
ax[0].hist(errs, bins=100)
ax[1].hist(errs_generic, bins=100)

popts = np.array(popts)

fig, ax = plt.subplots(2, 2)
ax[0, 0].plot(p_space, popts[:, 0])
ax[0, 1].plot(p_space, popts[:, 1])
ax[1, 0].plot(p_space, popts[:, 2])
ax[1, 1].plot(p_space, popts[:, 3])

plt.show()
