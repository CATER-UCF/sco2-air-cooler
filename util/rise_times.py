import pandas as pd

file_name = './data/time_series_step_changes_p8_e7.csv'
df = pd.read_csv(file_name)


def calculate_times(var_name, rise_pct, step_down_start, step_up_start):

    var = df[var_name]
    time = df['time']

    idx_down = df.loc[time == step_down_start].index[0]
    idx_up = df.loc[time == step_up_start].index[0]

    start_var_down = var[idx_down - 1]
    start_var_up = var[idx_up - 1]

    target_rise = rise_pct * (max(var) - start_var_up) + start_var_up
    target_fall = rise_pct * (min(var) - start_var_down) + start_var_down

    df_rise = df.loc[var > target_rise]
    df_fall = df.loc[var < target_fall]

    return df_rise['time'].iloc[0] - step_up_start, df_fall['time'].iloc[0] - step_down_start


temp_rise, temp_fall = calculate_times('temperature_tube_out_55', 0.9, 300, 900)
dense_rise, dense_fall = calculate_times('density_tube_out_55', 0.9, 900, 300)

print(f'Temperature rise time (90%): {temp_rise}   Fall time: {temp_fall}')
print(f'Density rise time (90%): {dense_rise}   Fall time: {dense_fall}')


temp_rise, temp_fall = calculate_times('temperature_tube_out_55', 0.95, 300, 900)
dense_rise, dense_fall = calculate_times('density_tube_out_55', 0.95, 900, 300)

print(f'Temperature rise time (95%): {temp_rise}   Fall time: {temp_fall}')
print(f'Density rise time (95%): {dense_rise}   Fall time: {dense_fall}')
