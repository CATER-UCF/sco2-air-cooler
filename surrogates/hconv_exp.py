"""
Surrogate model using an exponential function of the form:

hconv = a + sum(b_i * exp(-c_i(T - mu) ** e_i))

...for each pressure. This seems to give a better fit. Though interpolating
between pressures might pose some difficulty.
"""
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from scipy import optimize
import matplotlib.pyplot as plt

npr = 30
nt = 99

df_star = pd.read_csv(f'./data/p{npr}_t_star.csv')
p_space = np.array(df_star['pressure']).reshape(-1, 1)
t_star = np.array(df_star['T*'])
poly3 = PolynomialFeatures(degree=3)

p3 = poly3.fit_transform(p_space)
t_star_model = LinearRegression().fit(p3, t_star)

fig, ax = plt.subplots(1)
ax.plot(p_space, t_star_model.predict(p3), label='Fit')
ax.plot(p_space, t_star, '.', label='Act.')
ax.set_xlabel('Pressure (Pa)')
ax.set_xlabel('T* (K)')
ax.legend()


def t_star_pred(p):
    p3 = poly3.fit_transform(np.array(p).reshape(-1, 1))
    return t_star_model.predict(p3)

df = pd.read_csv(f'./data/DOE_p{npr}_t{nt}.csv')
df = df[df['temperature'] >= 300]

fig, ax = plt.subplots(2)

for i in range(npr):
    pi = p_space[i]
    df_tmp = df[abs(df['pressure'] - pi) < 0.01].sort_values('temperature')
    ti = np.array(df_tmp['temperature'])
    hi = np.array(df_tmp['hconv'])
    dpi = np.array(df_tmp['dP_over_l'])

    ax[0].plot(ti, hi, label=f'p={pi}')


def exp_kernel(x, a, mu, b1, c1, b2, c2, b3, c3):
    return a + b1 * np.exp(-c1 * np.abs(x - mu) ** 1) \
           + b2 * np.exp(-c2 * np.abs(x - mu) ** 0.5) \
           + b3 * np.exp(-c3 * np.abs(x - mu) ** 2)


popts = []
plabels = ['a', 'mu', 'b1', 'c1', 'b2', 'c2', 'b3', 'c3']

for i in range(npr):

    pi = p_space[i]
    df_tmp = df[abs(df['pressure'] - pi) < 0.01].sort_values('temperature')
    ti = np.array(df_tmp['temperature'])
    hi = np.array(df_tmp['hconv'])

    a = min(hi)
    b = (max(hi) - min(hi)) / 3
    mu = t_star_pred(p_space[i])[0]

    p_init = (a, mu, b, 1, b, 1, b, 1)

    try:
        popt, pcov = optimize.curve_fit(exp_kernel, ti, hi, p0=p_init)
    except RuntimeError:
        popt = p_init
        print('Optimal parameters not found!!!')

    popts.append(popt)
    hfit = exp_kernel(ti, *popt)
    hinit = exp_kernel(ti, *p_init)

    ax[1].plot(ti, hfit, label=f'hconv fit, P={pi}')

    if i % 10 == 0:
        fig, ax1 = plt.subplots(1)
        ax1.plot(ti, hi, label='hconv')
        ax1.plot(ti, hfit, label='hconv fit')
        ax1.plot(ti, hinit, label='hconv init')
        ax1.legend()

p_arr = np.array(popts)
fig, ax = plt.subplots(3, 3)

for i in range(p_arr.shape[1]):
    ax[i % 3, i // 3].plot(p_space, p_arr[:, i])
    ax[i % 3, i // 3].set_title(plabels[i])

plt.show()
