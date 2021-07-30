"""
Surrogate model of hconv using linear regression. This doesn't work very well...
"""
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt

npr = 30
nt = 99

# First, we want a model that tells us if which side of the ridge we're on
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

fig, ax = plt.subplots(2)

for i in range(npr):
    pi = p_space[i]
    df_tmp = df[abs(df['pressure'] - pi) < 0.01].sort_values('temperature')
    ti = np.array(df_tmp['temperature'])
    hi = np.array(df_tmp['hconv'])
    dpi = np.array(df_tmp['dP_over_l'])

    ax[0].plot(ti, hi, label=f'p={pi}')
    ax[1].plot(ti, dpi, label=f'p={pi}')


df['side'] = df.apply(lambda x: t_star_pred(x.pressure) > x.temperature, axis=1)
df['tst_inv'] = df.apply(lambda x: 1 / (x.temperature - t_star_pred(x.pressure) + 0.001), axis=1)

df_gr = df[df['side'] == False]
df_gr = df_gr[df_gr['temperature'] < 310]


print(df_gr['tst_inv'])

X = np.zeros((len(df_gr), 3))

temps = df_gr['temperature']
X[:, 0] = df_gr['temperature']
X[:, 1] = df_gr['pressure']
X[:, 2] = df_gr['tst_inv']

poly7 = PolynomialFeatures(degree=5)
XX = poly7.fit_transform(X)

model = LinearRegression().fit(XX, df_gr['hconv'])
hconv_gr = model.predict(XX)
print(model.coef_)
print(model.intercept_)

fig, ax = plt.subplots(2)

ax[0].plot(temps, df_gr['hconv'], '.', label='Act')
ax[0].plot(temps, hconv_gr, '.', label='Fit')
ax[1].plot(temps, df_gr['tst_inv'], '.')

ax[0].legend()
ax[1].legend()
plt.show()
