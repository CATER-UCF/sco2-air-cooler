"""
Flowsheet for generating a surrogate model of hconv and dP as functions of
fluid state.

Work in progress...
"""
import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pe
import pyomo as pyo
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties import swco2
from models import TubeSurrogate
import pandas as pd

N_DOE = 30


def make_model(sco2_t, sco2_p, sco2_mol_flow, dyn=True, n_pts=10):

    m = pe.ConcreteModel()
    if dyn:
        m.fs = FlowsheetBlock(default={"dynamic": True,
                                       "time_set": [0, n_pts],
                                       "time_units": pe.units.s})
    else:
        m.fs = FlowsheetBlock()

    m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
    m.fs.feed = TubeSurrogate(default={"property_package": m.fs.prop_sco2})

    m.fs.feed.add_geometry()
    m.fs.feed.use_gnielinski()
    m.fs.feed.use_churchill()

    if dyn:
        m.discretizer = pe.TransformationFactory('dae.finite_difference')
        m.discretizer.apply_to(m, nfe=n_pts - 1, wrt=m.fs.time, scheme="BACKWARD")

    m.fs.feed.tube_inner_perimeter = 0.0275 * 3.14159
    m.fs.feed.tube_length = 195
    m.fs.feed.number_of_tubes = 1160
    m.fs.feed.tube_hydraulic_diameter = 0.0275
    m.fs.feed.tube_flow_area = 3.14159 * (0.0275 / 2) ** 2

    sco2_enthalpy = swco2.htpx(T=sco2_t * pe.units.K, P=sco2_p * pe.units.Pa)

    m.fs.feed.outlet.flow_mol[:].fix(sco2_mol_flow)
    m.fs.feed.outlet.pressure[:].fix(sco2_p)
    m.fs.feed.outlet.enth_mol[:].fix(sco2_enthalpy)

    return m


# Full-factorial DOE over a range of temperatures and pressures
t_min = 300
t_max = 390
p_min = 7450000
p_max = 7760000
p_space = np.linspace(p_min, p_max, N_DOE)
t_star = np.zeros_like(p_space)


solver = pe.SolverFactory('ipopt')
solver.options = {
            "tol": 1e-6,
            "linear_solver": "ma27",
            "max_iter": 500,
        }

"""
For each pressure, we want to find the temperature which maximizes hconv (T*). 
So setup an optimization problem for each of these.
"""
def get_hconv(m):
    return m.fs.feed.hconv_tube[0]


def t_lower_bound(m):
    return m.fs.feed.properties[0].temperature >= t_min


def t_upper_bound(m):
    return m.fs.feed.properties[0].temperature <= t_max


for i, p in enumerate(p_space):

    # Setup and solve the model
    m = make_model(384.35, p, 13896.84163, dyn=False)
    m.fs.feed.initialize()
    solver.solve(m)

    # Find the critical temperature for this pressure
    m.fs.feed.outlet.enth_mol[0].unfix()
    m.obj = pe.Objective(rule=get_hconv, sense=pe.maximize)
    m.lb = pe.Constraint(rule=t_lower_bound)
    m.ub = pe.Constraint(rule=t_upper_bound)

    solver.solve(m)
    t_star[i] = pe.value(m.fs.feed.properties[0].temperature)
    print(f'T* found... for P={p}')

"""
Now setup our DOE and make sure each T*, P pair is included.
"""
cheb1 = np.polynomial.chebyshev.chebpts1(N_DOE // 2)
cheb1 = cheb1 / cheb1[-1]

cheb2 = np.polynomial.chebyshev.chebpts1(N_DOE // 2 + 1 + N_DOE % 2)
cheb2 = cheb2 / cheb2[-1]

presss = np.zeros(N_DOE ** 2)
temps = np.zeros(N_DOE ** 2)

for i, p in enumerate(p_space):
    t_range_1 = t_min + (cheb1 + 1) / 2 * (t_star[i] - t_min)
    t_range_2 = t_star[i] + (cheb2 + 1) / 2 * (t_max - t_star[i])
    t_full = np.concatenate((t_range_1, t_range_2[1:]))
    temps[i * N_DOE: (i + 1) * N_DOE] = t_full
    presss[i * N_DOE: (i + 1) * N_DOE] = np.ones_like(t_full) * p

m = make_model(384.35, 7751362.5, 13896.84163, dyn=True, n_pts=N_DOE ** 2)

for i, t in enumerate(m.fs.time):
    sco2_enthalpy = swco2.htpx(T=temps[i] * pe.units.K, P=presss[i] * pe.units.Pa)
    m.fs.feed.outlet.pressure[t].fix(presss[i])
    m.fs.feed.outlet.enth_mol[t].fix(sco2_enthalpy)

m.fs.feed.initialize()
solver = pe.SolverFactory('ipopt')

solver.options = {
            "tol": 1e-6,
            "linear_solver": "ma27",
            "max_iter": 500,
        }

solver.solve(m, tee=True)

cp_mol = np.array(pe.value(m.fs.feed.properties[:].cp_mol))
hconv = np.array(pe.value(m.fs.feed.hconv_tube[:]))
HTC = np.reshape(hconv, (N_DOE, N_DOE))

"""
fig, ax = plt.subplots(1)
cntr = ax.contour(temps, presss, HTC)
cb = plt.colorbar(cntr)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (Pa)')
ax.set_title('HTC as a Function of T and P')
"""

fig, ax = plt.subplots(2)
ax[0].plot(temps, hconv, '.')
ax[1].plot(presss, hconv, '.')
ax[0].set_xlabel('Temperature (K)')
ax[1].set_xlabel('Pressure (Pa)')


dp = np.array(pe.value(m.fs.feed.dP_over_l[:])) * 195
DP = np.reshape(dp, (N_DOE, N_DOE))

"""
fig, ax = plt.subplots(1)
cntr = ax.contour(temps, presss, DP)
cb = plt.colorbar(cntr)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (Pa)')
ax.set_title('Pressure Loss as a Function of T and P')
"""

fig, ax = plt.subplots(2)
ax[0].plot(temps, dp, '.')
ax[1].plot(presss, dp, '.')
ax[0].set_xlabel('Temperature (K)')
ax[1].set_xlabel('Pressure (Pa)')


df = pd.DataFrame(data={
    'temperature': temps,
    'pressure': presss,
    'hconv': hconv,
    'dP_over_l': dp / 195,
    'cp_mol': cp_mol
})
df.to_csv('./data/DOE_30.csv')

plt.show()
