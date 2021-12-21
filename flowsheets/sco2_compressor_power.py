import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pe
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties import swco2
from idaes.generic_models.unit_models import Compressor
import pandas as pd


df1 = pd.read_csv('./data/steady_state_vs_temperature_low.csv')
df2 = pd.read_csv('./data/steady_state_vs_temperature_crit.csv')
df3 = pd.read_csv('./data/steady_state_vs_temperature_high.csv')
frames = [df1, df2, df3]

df = pd.concat(frames, ignore_index=True)

df = pd.read_csv('./data/combined.csv')

temperatures = df['temperature_tube_out_55']
pressures = df['pressure_tube_out_55']

n_pts = len(temperatures)
m = pe.ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": True,
                               "time_set": [0, 1],
                               "time_units": pe.units.s})

m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
m.fs.compressor = Compressor(default={"property_package": m.fs.prop_sco2,
                                      "dynamic": False})

m.discretizer = pe.TransformationFactory('dae.finite_difference')
m.discretizer.apply_to(m, nfe=n_pts - 1, wrt=m.fs.time, scheme="BACKWARD")

# Assume 20MPa and 70% efficiency
m.fs.compressor.outlet.pressure[:].fix(20e6)
m.fs.compressor.efficiency_isentropic.fix(0.7)

# Inlet conditions
m.fs.compressor.inlet.flow_mol[:].fix(13896.84163)
for i, t in enumerate(m.fs.time):
    enth = swco2.htpx(T=temperatures[i] * pe.units.K, P=pressures[i] * pe.units.Pa)
    m.fs.compressor.inlet.pressure[t].fix(pressures[i])
    m.fs.compressor.inlet.enth_mol[t].fix(enth)

m.fs.compressor.initialize()

solver = pe.SolverFactory('ipopt')
solver.options = {
    "tol": 1e-6,
    "linear_solver": "ma27",
    "max_iter": 500,
}

solver.solve(m, tee=True)

work = pe.value(m.fs.compressor.work_mechanical[:])
df['compressor_power'] = np.array(work) * 1e-6

df.to_csv('./data/combined.csv', index=None)
