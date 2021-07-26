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


m = make_model(384.35, 7751362.5, 13896.84163, dyn=True, n_pts=N_DOE ** 2)

# Full-factorial DOE over a range of temperatures and pressures
t_min = 300
t_max = 310
p_min = 7450000
p_max = 7760000
cheb = np.polynomial.chebyshev.chebpts1(N_DOE)
cheb = cheb / cheb[-1]

t = t_min + (cheb + 1) / 2 * (t_max - t_min)
p = p_min + (cheb + 1) / 2 * (p_max - p_min)
t_sco2, p_sco2 = np.meshgrid(t, p)

t_flat = t_sco2.flatten()
p_flat = p_sco2.flatten()

for i, t in enumerate(m.fs.time):
    sco2_enthalpy = swco2.htpx(T=t_flat[i] * pe.units.K, P=p_flat[i] * pe.units.Pa)
    m.fs.feed.outlet.pressure[t].fix(p_flat[i])
    m.fs.feed.outlet.enth_mol[t].fix(sco2_enthalpy)

m.fs.feed.initialize()
solver = pe.SolverFactory('ipopt')

solver.options = {
            "tol": 1e-6,
            "linear_solver": "ma27",
            "max_iter": 500,
        }

solver.solve(m, tee=True)

hconv = np.array(pe.value(m.fs.feed.hconv_tube[:]))
HTC = np.reshape(hconv, (N_DOE, N_DOE))

fig, ax = plt.subplots(1)
cntr = ax.contourf(t_sco2, p_sco2, HTC)
cb = plt.colorbar(cntr)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (Pa)')
ax.set_title('HTC as a Function of T and P')

fig, ax = plt.subplots(2)
ax[0].plot(t_flat, hconv, '.')
ax[1].plot(p_flat, hconv, '.')
ax[0].set_xlabel('Temperature (K)')
ax[1].set_xlabel('Pressure (Pa)')


dp = np.array(pe.value(m.fs.feed.dP_over_l[:])) * 195
DP = np.reshape(dp, (N_DOE, N_DOE))

fig, ax = plt.subplots(1)
cntr = ax.contourf(t_sco2, p_sco2, DP)
cb = plt.colorbar(cntr)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (Pa)')
ax.set_title('Pressure Loss as a Function of T and P')

fig, ax = plt.subplots(2)
ax[0].plot(t_flat, dp, '.')
ax[1].plot(p_flat, dp, '.')
ax[0].set_xlabel('Temperature (K)')
ax[1].set_xlabel('Pressure (Pa)')


df = pd.DataFrame(data={
    'temperature': t_flat,
    'pressure': p_flat,
    'hconv': hconv,
    'dP_over_l': dp / 195
})
df.to_csv('./data/DOE_30.csv')

plt.show()
