"""
Surrogate model using an exponential function of the form:

hconv = a + sum(b_i * exp(-c_i(T - mu) ** e_i))

...at two different pressures, then linearly interpolating between them. This
produces a poor fit since T* is different for each pressure. It might be
possible to improve this if we interpolate between (P - P*) somehow.
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


def exp_kernel1(x, a, mu, b1, c1, b2, c2, b3, c3):
    return a + b1 * np.exp(-c1 * np.abs(x - mu) ** 1) \
           + b2 * np.exp(-c2 * np.abs(x - mu) ** 1.5) \
           + b3 * np.exp(-c3 * np.abs(x - mu) ** 2)


def exp_kernel2(x, a, mu, b1, c1, b2, c2, b3, c3):
    return a + b1 * np.exp(-c1 * np.abs(x - mu) ** 4) \
           + b2 * np.exp(-c2 * np.abs(x - mu) ** 3) \
           + b3 * np.exp(-c3 * np.abs(x - mu) ** 2)


def make_fit(pi, kernel, title):

    df_tmp = df[abs(df['pressure'] - pi) < 0.01].sort_values('temperature')
    ti = np.array(df_tmp['temperature'])
    hi = np.array(df_tmp['hconv'])

    a = min(hi)
    b = (max(hi) - min(hi)) / 3
    mu = t_star_pred(pi)[0]

    p_init = (a, mu, b, 1, b, 1, b, 1)
    popt, pcov = optimize.curve_fit(kernel, ti, hi, p0=p_init)

    hfit = kernel(ti, *popt)
    hinit = kernel(ti, *p_init)

    fig, ax1 = plt.subplots(2)
    ax1[0].plot(ti, hi, label='hconv')
    ax1[0].plot(ti, hfit, label='hconv fit')
    ax1[0].plot(ti, hinit, label='hconv init')
    ax1[0].legend()
    ax1[0].set_title(title)
    ax1[1].plot(ti, hfit - hi, label='Error')
    ax1[1].legend()

    return popt


p1 = min(p_space)[0]
p2 = max(p_space)[0]
popt1 = make_fit(p1, exp_kernel1, 'kernel 1')
popt2 = make_fit(p2, exp_kernel2, 'kernel 2')

errors = []
pct_errors = []

for i in range(npr):

    pi = p_space[i][0]
    df_tmp = df[abs(df['pressure'] - pi) < 0.01].sort_values('temperature')
    ti = np.array(df_tmp['temperature'])
    hi = np.array(df_tmp['hconv'])

    exp1 = exp_kernel1(ti, *popt1)
    exp2 = exp_kernel2(ti, *popt2)

    hfit = exp1 + (exp2 - exp1) * (pi - p1) / (p2 - p1)

    errors.append(hfit - hi)
    pct_errors.append(100 * (hfit - hi) / hi)

    if abs(pi - p1) < 0.01 or abs(pi - p2) < 0.01:
        ax[1].plot(ti, hfit, label=f'hconv fit, P={pi}')


errs = np.concatenate(errors, axis=None)
pct_errs = np.concatenate(pct_errors, axis=None)
print(f'mean: {np.mean(errs)}, std: {np.std(errs)}, min: {np.min(errs)}, max: {np.max(errs)}')
print(f'mean: {np.mean(pct_errs)}, std: {np.std(pct_errs)}, min: {np.min(pct_errs)}, max: {np.max(pct_errs)}')

fig, ax = plt.subplots(2)
ax[0].hist(errs, bins=30)
ax[1].hist(pct_errs, bins=30)

plt.show()
