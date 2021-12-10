import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

matplotlib.rcParams.update({'font.size': 12})


def plot_response(data_file, n_passes=8, n_elements=7,
                  show=False, image_file=None,
                  loc0='upper left', loc1='upper left',
                  response_var='temperature_tube_out', response_label='CO2 Out',
                  response_axis_label='Temperature ($^\circ$C)'):

    # Read in the results data
    df = pd.read_csv(data_file)
    n_eles = n_passes * n_elements

    time = df['time']
    air_temp = df[f'temperature_shell_in_{n_eles - 1}']
    response = df[f'{response_var}_{n_eles - 1}']
    if 'temperature' in response_var:
        response -= 273.15

    fig = plt.figure(constrained_layout=True, figsize=(6, 5))
    gs = fig.add_gridspec(2, 1)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])

    for ax in [ax0, ax1]:
        ax.set_axisbelow(True)
        ax.xaxis.grid(color='gray', linestyle='dashed')
        ax.yaxis.grid(color='gray', linestyle='dashed')

    ax0.plot(time, response, 'k', ls='solid', label=response_label)
    ax1.plot(time, air_temp - 273.15, 'k', ls='dashdot', label='Air In')

    #ax0.get_xaxis().set_visible(False)
    ax0.set_ylabel(response_axis_label)

    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('Temperature ($^\circ$C)')

    ax0.legend(loc=loc0)
    ax1.legend(loc=loc1)

    if show:
        plt.show()

    if image_file is not None:
        fig.savefig(image_file)


plot_response('./data/time_series_step_changes_p8_e7.csv',
              image_file='./images/temperature_response_step_changes.png')

plot_response('./data/time_series_ramp_down_p8_e7.csv',
              loc0='upper right', loc1='upper right',
              image_file='./images/temperature_response_ramp_down.png')

plot_response('./data/time_series_ramp_up_p8_e7.csv',
              loc0='upper left', loc1='upper left',
              image_file='./images/temperature_response_ramp_up.png')

plot_response('./data/time_series_step_changes_p8_e7.csv',
              loc0='upper right',
              response_var='density_tube_out',
              response_axis_label='Density (kg/m$^3$)',
              image_file='./images/density_response_step_changes.png')

plot_response('./data/time_series_ramp_down_p8_e7.csv',
              loc0='upper left', loc1='upper right',
              response_var='density_tube_out',
              response_axis_label='Density (kg/m$^3$)',
              image_file='./images/density_response_ramp_down.png')

plot_response('./data/time_series_ramp_up_p8_e7.csv',
              loc0='upper right', loc1='upper left',
              response_var='density_tube_out',
              response_axis_label='Density (kg/m$^3$)',
              image_file='./images/density_response_ramp_up.png')
