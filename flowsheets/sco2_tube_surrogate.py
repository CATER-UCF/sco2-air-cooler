"""
Flowsheet for generating a surrogate model of hconv and dP as functions of
fluid state.

Work in progress...
"""
import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pe
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties import swco2
from models import TubeSurrogate

N_DOE = 20


def make_model(sco2_t, sco2_p, sco2_mol_flow):

    m = pe.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   "time_set": [0, N_DOE ** 2],
                                   "time_units": pe.units.s})

    m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
    m.fs.feed = TubeSurrogate(default={"property_package": m.fs.prop_sco2})

    m.fs.feed.add_geometry()
    m.fs.feed.add_hconv_eqs()

    m.discretizer = pe.TransformationFactory('dae.finite_difference')
    m.discretizer.apply_to(m, nfe=N_DOE ** 2 - 1, wrt=m.fs.time, scheme="BACKWARD")

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


m = make_model(384.35, 7751362.5, 13896.84163)

m.fs.feed.initialize()
solver = pe.SolverFactory('ipopt')


# Full-factorial DOE over a range of temperatures and pressures
t_min = 300
t_max = 400
p_min = 7450000
p_max = 7760000
t = np.linspace(t_min, t_max, N_DOE)
p = np.linspace(p_min, p_max, N_DOE)
t_sco2, p_sco2 = np.meshgrid(t, p)

t_flat = t_sco2.flatten()
p_flat = p_sco2.flatten()

for i, t in enumerate(m.fs.time):
    sco2_enthalpy = swco2.htpx(T=t_flat[i] * pe.units.K, P=p_flat[i] * pe.units.Pa)
    m.fs.feed.outlet.pressure[t].fix(p_flat[i])
    m.fs.feed.outlet.enth_mol[t].fix(sco2_enthalpy)

solver.options = {
            "tol": 1e-6,
            "linear_solver": "ma27",
            "max_iter": 500,
        }

solver.solve(m, tee=True)

h = pe.value(m.fs.feed.hconv_tube[:])
HTC = np.reshape(h, (N_DOE, N_DOE))
plt.contourf(t_sco2, p_sco2, HTC)
cb = plt.colorbar()
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
plt.title('HTC as a Function of T and P')
plt.show()
