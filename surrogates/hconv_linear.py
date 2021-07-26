import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt

df = pd.read_csv('./data/DOE_30.csv')

X = np.zeros((len(df), 2))
X[:, 0] = df['temperature']
X[:, 1] = df['pressure']

poly = PolynomialFeatures(degree=3)
XX = poly.fit_transform(X)

model = LinearRegression().fit(XX, df['hconv'])
hconv_fit = model.predict(XX)
print(model.coef_)
print(model.intercept_)

fig, ax = plt.subplots(2)
ax[0].plot(df['temperature'], df['hconv'], '.', label='Act')
ax[1].plot(df['pressure'], df['hconv'], '.', label='Act')
ax[0].plot(df['temperature'], hconv_fit, '.', label='Fit')
ax[1].plot(df['pressure'], hconv_fit, '.', label='Fit')

ax[0].legend()
ax[1].legend()


fig, ax = plt.subplots(2)
temps = np.array(df['temperature'])
hconvs = np.array(df['hconv'])
cps = np.array(df['cp_mol'])
presss = np.array(df['pressure'])
for i in range(30):
    ti = temps[i * 30: i * 30 + 30]
    hi = hconvs[i * 30: i * 30 + 30]
    cpi = cps[i * 30: i * 30 + 30]
    pi = presss[i * 30]
    ax[0].plot(ti, hi, label=f'p={pi}')
    ax[1].plot(ti, cpi, label=f'p={pi}')

#ax.legend()
plt.show()
