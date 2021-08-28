import pandas as pd
import pyomo.environ as pe
import numpy as np


def get_arr(expr):
    return np.array(pe.value(expr))


def write_csv(file_name, elements):
    data = {'time': np.array([t for t in elements[0].flowsheet().time])}
    for idx, e in enumerate(elements):
        t_tube = get_arr(e.tube.properties_in[:].temperature)
        p_tube = get_arr(e.tube.properties_in[:].pressure)
        t_shell = get_arr(e.shell.properties_in[:].temperature)
        heat_duty = get_arr(e.heat_duty[:]) * 1e-6
        data[f'temperature_tube_{idx}'] = t_tube
        data[f'pressure_tube_{idx}'] = p_tube
        data[f'temperature_shell_{idx}'] = t_shell
        data[f'heat_duty_{idx + 1}'] = heat_duty

    idx = len(elements)
    data[f'temperature_tube_{idx}'] = get_arr(elements[-1].tube
                                              .properties_out[:].temperature)
    data[f'pressure_tube_{idx}'] = get_arr(elements[-1].tube
                                           .properties_out[:].pressure)
    data[f'temperature_shell_{idx}'] = get_arr(elements[-1].shell
                                               .properties_out[:].temperature)

    df = pd.DataFrame(data=data)
    df.to_csv(file_name, index=None)
