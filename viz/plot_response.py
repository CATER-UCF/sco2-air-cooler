import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

matplotlib.rcParams.update({'font.size': 14})


def plot_2d_steady_state(data_file, n_passes=8, n_elements=7,
                         show=False, image_file=None):

    # Read in the results data
    df = pd.read_csv(data_file)
    n_eles = n_passes * n_elements

    time = df['time']
    air_temp = df[f'temperature_shell_{n_eles - 1}']
    sco2_temp = df[f'temperature_tube_{n_eles}']

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(time, sco2_temp - 273.15, 'b', label='CO2 Out')
    ax[1].plot(time, air_temp - 273.15, 'r', label='Air In')

    [axs.legend(loc='upper left') for axs in ax]
    ax[1].set_xlabel('time (s)')
    ax[0].get_xaxis().set_visible(False)

    fig.suptitle('Temperature Response (C)')

    if show:
        plt.show()

    if image_file is not None:
        fig.savefig(image_file)


plot_2d_steady_state('./data/time_series_p8_e7.csv',
                     n_passes=7, n_elements=8, show=False,
                     image_file='./images/temp_response_p8_e7.png')
