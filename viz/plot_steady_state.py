import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
import pandas as pd
import numpy as np

data_file = './data/steady_state_vs_temperature_p8_e7.csv'

df = pd.read_csv(data_file)

# sCO2 exit temperature, density vs air temperature
air_temp = df['temperature_shell_in_55'] - 273.15
sco2_temp = df['temperature_tube_out_55'] - 273.15
sco2_press = df['pressure_tube_out_55'] * 1e-5
sco2_rho = df['density_tube_out_55']
sco2_hconv = df['hconv_tube_55']

matplotlib.rcParams.update({'font.size': 12})
fig = plt.figure(constrained_layout=True, figsize=(6, 5))
gs = fig.add_gridspec(2, 1)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])

for ax in [ax0, ax1]:
    ax.set_axisbelow(True)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(color='gray', linestyle='dashed')

ax0.plot(air_temp, sco2_temp, 'k')
ax1.plot(air_temp, sco2_rho, 'k')

ax0.set_ylabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Density (kg/m$^3$)')
ax1.set_xlabel('Air Temperature ($^\circ$C)')

fig.savefig('./images/steady_state_vs_temperature.png')


# ... vs tube length
rc('font', **{'family': 'sans-serif', 'sans-serif': ['DejaVu Sans'], 'size': 10})
rc('mathtext', **{'default': 'regular'})
fig, ax = plt.subplots(2, 2, figsize=(8, 6))

for i in range(2):
    for j in range(2):
        ax[i][j].set_axisbelow(True)
        ax[i][j].xaxis.grid(color='gray', linestyle='dashed')
        ax[i][j].yaxis.grid(color='gray', linestyle='dashed')

temp = np.empty((56, 1))
press = np.empty((56, 1))
dens = np.empty((56, 1))
hconv = np.empty((56, 1))
tube_length = 195
x = np.linspace(0, tube_length, 56)

d = {'air_temp': [273.15, 288.15, 303.15],
     'ls': ['solid', 'dashed', 'dotted'],
     'color': ['k', '0.25', '0.10']}

for k, a in enumerate(d['air_temp']):
    idx = df.loc[abs(df['temperature_shell_in_55'] - a) < 0.01].index[0]
    for j in range(56):
        temp[j] = df.iloc[idx][f'temperature_tube_out_{j}'] - 273.15
        press[j] = df.iloc[idx][f'pressure_tube_out_{j}'] * 1e-5
        dens[j] = df.iloc[idx][f'density_tube_out_{j}']
        hconv[j] = df.iloc[idx][f'hconv_tube_{j}']

    c = d['color'][k]
    ls = d['ls'][k]
    ax[0][0].plot(x, temp, c, ls=ls, label=f'$T_{{air}}$: {int(a - 273.15)}$^\circ$C')
    ax[0][1].plot(x, press, c, ls=ls)
    ax[1][0].plot(x, dens, c, ls=ls)
    ax[1][1].plot(x, hconv, c, ls=ls)

ax[0][0].legend(loc='upper right')
ax[1][0].set_xlabel('Tube Position (m)')
ax[1][1].set_xlabel('Tube Position (m)')

ax[0][0].set_ylabel('Temperature ($^\circ$C)')
ax[0][1].set_ylabel('Pressure (Bar)')
ax[1][0].set_ylabel('Density (W/m$^2$K)')
ax[1][1].set_ylabel('hconv (W/m$^2$K)')

plt.subplots_adjust(wspace=0.4, top=0.95, right=0.95, left=0.1)
fig.savefig('./images/temperature_profiles.png')
plt.show()
