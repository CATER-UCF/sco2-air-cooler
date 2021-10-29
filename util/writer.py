import pandas as pd
import pyomo.environ as pe
import numpy as np


def get_arr(expr):
    return np.array(pe.value(expr))


def write_csv(file_name, elements):

    data = {'time': np.array([t for t in elements[0].flowsheet().time])}
    for idx, e in enumerate(elements):

        t_tube_in = get_arr(e.tube.properties_in[:].temperature)
        t_tube_out = get_arr(e.tube.properties_out[:].temperature)
        p_tube_in = get_arr(e.tube.properties_in[:].pressure)
        p_tube_out = get_arr(e.tube.properties_out[:].pressure)
        rho_tube_in = get_arr(e.tube.properties_in[:].dens_mass)
        rho_tube_out = get_arr(e.tube.properties_out[:].dens_mass)
        t_shell_in = get_arr(e.shell.properties_in[:].temperature)
        t_shell_out = get_arr(e.shell.properties_out[:].temperature)
        heat_duty = get_arr(e.heat_duty[:]) * 1e-6

        data[f'temperature_tube_in_{idx}'] = t_tube_in
        data[f'temperature_tube_out_{idx}'] = t_tube_out
        data[f'pressure_tube_in_{idx}'] = p_tube_in
        data[f'pressure_tube_out_{idx}'] = p_tube_out
        data[f'density_tube_in_{idx}'] = rho_tube_in
        data[f'density_tube_out_{idx}'] = rho_tube_out
        data[f'temperature_shell_in_{idx}'] = t_shell_in
        data[f'temperature_shell_out_{idx}'] = t_shell_out
        data[f'heat_duty_{idx}'] = heat_duty

    df = pd.DataFrame(data=data)
    df.to_csv(file_name, index=None)
