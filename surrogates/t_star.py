"""
Polynomial fit for critical temperature as a function of pressure.

T* = 2.55996364e+02 + 7.15290022e-06 * P -8.52339353e-14 * P ** 2
"""
import pandas as pd
import numpy as np

npr = 30

df_star = pd.read_csv(f'./data/p{npr}_t_star.csv')
p = np.array(df_star['pressure'])
t_star = np.array(df_star['T*'])

coeffs = np.polyfit(p, t_star, 2)

print(coeffs)
