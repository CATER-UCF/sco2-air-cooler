import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np


matplotlib.rcParams.update({'font.size': 12})
data_file = './data/steady_state_vs_temperature_p8_e7.csv'

df = pd.read_csv(data_file)

# sCO2 exit temperature, density vs air temperature
air_temp = df['temperature_shell_in_55'] - 273.15
sco2_temp = df['temperature_tube_out_55'] - 273.15
sco2_rho = df['density_tube_out_55']

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

ax[0].plot(air_temp, sco2_temp, 'b', label='CO2 Out')
ax[1].plot(air_temp, sco2_rho, 'g', label='CO2 Out')


ax[0].set_ylabel('Temperature ($^\circ$C)')
ax[1].set_ylabel('Density (kg/m$^3$)')

ax[0].set_xlabel('Air Temperature ($^\circ$C)')
ax[1].set_xlabel('Air Temperature ($^\circ$C)')

ax[0].legend()
ax[1].legend()

fig.savefig('./images/steady_state_vs_temperature.png')


# ... vs tube length
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

temp = np.empty((55, 1))
dens = np.empty((55, 1))
tube_length = 195
x = np.linspace(0, tube_length, 55)

air_temps = [273.15, 288.15, 303.15]
for k, a in enumerate(air_temps):
    idx = df.loc[abs(df['temperature_shell_in_55'] - a) < 0.01].index[0]
    for j in range(55):
        temp[j] = df.iloc[idx][f'temperature_tube_in_{j}'] - 273.15
        dens[j] = df.iloc[idx][f'density_tube_in_{j}']

    ax[0].plot(x, temp, label=f'Air Temperature: {int(a - 273.15)}$^\circ$C')
    ax[1].plot(x, dens, label=f'Air Temperature: {int(a - 273.15)}$^\circ$C')

ax[0].legend()
ax[0].set_xlabel('Tube Position (m)')
ax[1].set_xlabel('Tube Position (m)')
ax[0].set_ylabel('Temperature ($^\circ$C)')
ax[1].set_ylabel('Density (kg/m$^3$)')

#ax[1].legend()

fig.savefig('./images/temperature_profiles.png')
