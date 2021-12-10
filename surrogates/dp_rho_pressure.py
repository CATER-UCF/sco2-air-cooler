"""
Surrogate model of dP / l using density:

dP / l =  8.82481388e+01 -1.81263728e+05 / rho +
         -2.98465078e-01 * rho + 1.92815541e-04 * rho ** 2

This yields an excellent fit with < 0.1% error. Adding a second property
state variable (i.e. pressure) makes for a marginal improvement.
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


def model(x, a0, a1, a2, a3):
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
    popt, _ = optimize.curve_fit(model, XX, Y, p0=p_init)
    return popt, XX, Y, temps


# "Generic" model valid for all pressures
popt_generic, _, _, _ = get_model(df)
print(popt_generic)

fig, ax = plt.subplots(2)
plt.subplots_adjust(hspace=0.5)

popts = []
errs = np.array([])
errs_generic = np.array([])

# Here, we try making a different model for each pressure
for p in p_space:
    df_tmp = df[abs(df['pressure'] - p[0]) < 0.01]
    popt, XX, Y, temps = get_model(df_tmp)
    Y_model = model(XX, *popt)
    Y_generic = model(XX, *popt_generic)
    err = (Y_model - Y) / Y
    err_generic = (Y_generic - Y) / Y
    ax[0].plot(temps, err, '.')
    ax[1].plot(temps, err_generic, '.')
    popts.append(popt)
    errs = np.append(errs, err)
    errs_generic = np.append(errs_generic, err_generic)

# Scatter plots...
ax[0].set_title('Pressure-Specific Models')
ax[1].set_title('Generic Model')
ax[0].set_ylabel('% Error')
ax[1].set_ylabel('% Error')
ax[0].set_xlabel('Temperature (K)')
ax[1].set_xlabel('Temperature (K)')

# Histograms...
fig, ax = plt.subplots(2)
plt.subplots_adjust(hspace=0.5)

ax[0].hist(errs, bins=100)
ax[1].hist(errs_generic, bins=100)
ax[0].set_title('Pressure-Specific Models')
ax[1].set_title('Generic Model')
ax[0].set_xlabel('% Error')
ax[1].set_xlabel('% Error')


popts = np.array(popts)

fig, ax = plt.subplots(2, 2)
ax[0, 0].plot(p_space, popts[:, 0])
ax[0, 1].plot(p_space, popts[:, 1])
ax[1, 0].plot(p_space, popts[:, 2])
ax[1, 1].plot(p_space, popts[:, 3])
fig.suptitle('Model Coefficients vs. Pressure')

plt.show()
