import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc


# File can be downloaded from: https://www.nrel.gov/gis/solar.html
df = pd.read_csv('./data/722050TYA_LADA.csv', skiprows=1)
df = df.rename(columns={'Date (MM/DD/YYYY)': 'Date', 'Time (HH:MM)': 'Time'})
temp_summer = df.loc[df['Date'] == '6/1/1996']['Dry-bulb (C)']
temp_winter = df.loc[df['Date'] == '1/1/1995']['Dry-bulb (C)']
times_summer = df.loc[df['Date'] == '6/1/1996']['Time']
times_winter = df.loc[df['Date'] == '1/1/1995']['Time']

# Times of day are in order, so this is ok...
print(times_summer, times_winter)
times = list(range(1, 25))

# Plot formatting
rc('font', **{'family': 'sans-serif', 'sans-serif': ['DejaVu Sans'], 'size': 8})
rc('mathtext', **{'default': 'regular'})

# Plot the temperatures
fig, ax = plt.subplots(constrained_layout=True, figsize=(3.255, 2))
ax.plot(times, temp_summer, 'k', linewidth=1, label='Summer')
ax.plot(times, temp_winter, 'k', linestyle='--', linewidth=1, label='Winter')

ticks = [3, 6, 9, 12, 15, 18, 21, 24]
labels = [str(t) + ':00' for t in ticks]
_ = ax.set_xticks(ticks, labels=labels, rotation=45)
ax.set_ylabel('Temperature ($^\circ$C)')
ax.legend(loc='upper left', handlelength=1)

fig.savefig('./images/ambient_temperatures.png', dpi=500)
plt.show()
