import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
import pandas as pd
import numpy as np


df = pd.read_csv('./data/combined.csv')

# sCO2 exit temperature, density vs air temperature
air_temp = df['temperature_shell_in_55'] - 273.15
sco2_temp = df['temperature_tube_out_55'] - 273.15
sco2_press = df['pressure_tube_out_55'] * 1e-5
sco2_rho = df['density_tube_out_55']
sco2_hconv = df['hconv_tube_55']
compressor_power = df['compressor_power']

matplotlib.rcParams.update({'font.size': 12})
fig = plt.figure(constrained_layout=True, figsize=(6, 7))
gs = fig.add_gridspec(3, 1)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[2, 0])

for ax in [ax0, ax1, ax2]:
    ax.set_axisbelow(True)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    ax.yaxis.grid(color='gray', linestyle='dashed')

ax0.plot(air_temp, sco2_temp, 'k')
ax1.plot(air_temp, sco2_rho, 'k')
ax2.plot(air_temp, compressor_power, 'k')

ax0.set_ylabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Density (kg/m$^3$)')
ax2.set_ylabel('Compressor Power (MW)')
ax2.set_xlabel('Air Temperature ($^\circ$C)')

fig.savefig('./images/steady_state_vs_temperature.png')


# ... vs tube length
rc('font', **{'family': 'sans-serif', 'sans-serif': ['DejaVu Sans'], 'size': 12})
rc('mathtext', **{'default': 'regular'})

fig1, ax1 = plt.subplots(2, 1, figsize=(6, 5))
fig2, ax2 = plt.subplots(2, 1, figsize=(6, 5))

for a in [ax1, ax2]:
    for i in range(2):
        a[i].set_axisbelow(True)
        a[i].xaxis.grid(color='gray', linestyle='dashed')
        a[i].yaxis.grid(color='gray', linestyle='dashed')

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
    ax1[0].plot(x, temp, c, ls=ls, label=f'$T_{{air}}$: {int(a - 273.15)}$^\circ$C')
    ax1[1].plot(x, press, c, ls=ls)
    ax2[0].plot(x, dens, c, ls=ls, label=f'$T_{{air}}$: {int(a - 273.15)}$^\circ$C')
    ax2[1].plot(x, hconv, c, ls=ls)

ax1[0].legend(loc='upper right')
ax2[0].legend(loc='upper left')

ax1[1].set_xlabel('Tube Position (m)')
ax2[1].set_xlabel('Tube Position (m)')

ax1[0].set_ylabel('Temperature ($^\circ$C)')
ax1[1].set_ylabel('Pressure (Bar)')
ax2[0].set_ylabel('Density (W/m$^2$K)')
ax2[1].set_ylabel('hconv (W/m$^2$K)')

fig1.subplots_adjust(wspace=0.4, top=0.95, right=0.95, left=0.15)
fig2.subplots_adjust(wspace=0.4, top=0.95, right=0.95, left=0.15)

fig1.savefig('./images/temperature_profiles1.png')
fig2.savefig('./images/temperature_profiles2.png')
plt.show()
