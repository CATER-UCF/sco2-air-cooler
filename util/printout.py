import pyomo.environ as pe
import pandas as pd


def print_results_0d(e, t=0):
    """
    Prints out the results from a 0D heat exchanger element.

    Parameters:
    -----------
        e : Heat exchanger unit model element
        t : time

    Returns:
    --------
        heat_duty : Heat transferred (MW)
    """
    t_tube_in = pe.value(e.tube.properties_in[t].temperature)
    p_tube_in = pe.value(e.tube.properties_in[t].pressure)
    h_tube_in = pe.value(e.tube.properties_in[t].enth_mol) / 18.01528

    t_tube_out = pe.value(e.tube.properties_out[t].temperature)
    p_tube_out = pe.value(e.tube.properties_out[t].pressure)
    h_tube_out = pe.value(e.tube.properties_out[t].enth_mol) / 18.01528

    t_shell_in = pe.value(e.shell.properties_in[t].temperature)
    p_shell_in = pe.value(e.shell.properties_in[t].pressure)

    t_shell_out = pe.value(e.shell.properties_out[t].temperature)
    p_shell_out = pe.value(e.shell.properties_out[t].pressure)

    heat_duty = pe.value(e.heat_duty[t])
    UA = pe.value(e.overall_heat_transfer_coefficient[t])

    crossflow_P = (t_tube_out - t_tube_in) / (t_shell_in - t_tube_in)
    crossflow_R = (t_shell_in - t_shell_out) / (t_tube_out - t_tube_in)

    print('--------Tube-----------'
          '')
    print('')
    print(f'T in: {t_tube_in}')
    print(f'T out: {t_tube_out}')
    print('')
    print(f'P in: {p_tube_in}')
    print(f'P out: {p_tube_out}')
    print('')
    print(f'H in: {h_tube_in}')
    print(f'H out: {h_tube_out}')
    print('')
    print('--------Shell-----------')
    print('')
    print(f'T in: {t_shell_in}')
    print(f'T out: {t_shell_out}')
    print('')
    print(f'P in: {p_shell_in}')
    print(f'P out: {p_shell_out}')
    print('')
    print(f'Overall HTC: {UA}')
    print(f'Heat duty (MW): {heat_duty * 1e-6}')
    print('')
    print(f'Crossflow R: {crossflow_R}')
    print(f'Crossflow P: {crossflow_P}')

    return heat_duty * 1e-6


def print_results_row(file_name, passes=8, elements=7, t=0):
    """
    Use for inspecting results from a transient simulation at one particular
    time point.

    Parameters:
    -----------
        file_name : .csv file name
        passes : number of tube passes
        elements : number of elements per pass
        t : time (row index)
    """
    df = pd.read_csv(file_name)

    last_element = passes * elements - 1

    temp_shell_out_cols = [f'temperature_shell_out_{i}' for i in range(elements)]
    heat_duty_cols = [f'heat_duty_{i+1}' for i in range(elements * passes)]

    temp_tube_in = df['temperature_tube_in_0'][t] - 273.15
    press_tube_in = df['pressure_tube_in_0'][t] * 1e-5
    temp_shell_in = df[f'temperature_shell_in_{last_element}'][t] - 273.15

    temp_tube_out = df[f'temperature_tube_out_{last_element}'][t] - 273.15
    press_tube_out = df[f'pressure_tube_out_{last_element}'][t] * 1e-5
    temp_shell_out = df[temp_shell_out_cols].mean(axis=1)[t] - 273.15
    heat_duty = df[heat_duty_cols].sum(axis=1)[t]

    print(f'CO2 Temperature In (C):    {temp_tube_in}')
    print(f'CO2 Pressure In (bar):     {press_tube_in}')
    print(f'Air Temperature In (C):    {temp_shell_in}')

    print(f'CO2 Temperature Out (C):   {temp_tube_out}')
    print(f'CO2 Pressure Out (bar):    {press_tube_out}')
    print(f'Air Temperature Out (C):   {temp_shell_out}')
    print(f'Heat Duty (MW):            {heat_duty}')


if __name__ == '__main__':
    print_results_row('./data/time_series_ramp_down_p8_e7.csv', t=0)
