"""
...
"""

import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pe
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import HeatExchanger, HeatExchangerFlowPattern, Mixer, Separator
from idaes.generic_models.unit_models.heat_exchanger import delta_temperature_lmtd_callback
from idaes.generic_models.properties import swco2
from idaes.power_generation.properties import FlueGasParameterBlock
import idaes.logger
from util import print_results_0d
from models import HeatExchangerElement
from flowsheets.code_gen import code_gen

logger = idaes.logger.getLogger('idaes')

logger.info('Creating model...')
m = pe.ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": True,
                               "time_set": [0, 300, 600, 900, 1200, 1500],
                               "time_units": pe.units.s})

m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
m.fs.prop_fluegas = FlueGasParameterBlock()

n_passes = 8
n_elements_per_pass = 7
n_elements = n_passes * n_elements_per_pass
m.fs.es = []

logger.info('Adding heat exchanger elements...')
for _ in range(n_elements):
    m.fs.es.append(HeatExchangerElement(default={
        "delta_temperature_callback": delta_temperature_lmtd_callback,
        "cold_side_name": "shell",
        "hot_side_name": "tube",
        "shell": {"property_package": m.fs.prop_fluegas,
                  "has_pressure_change": False},
        "tube": {"property_package": m.fs.prop_sco2,
                 "has_pressure_change": False},
        "flow_pattern": HeatExchangerFlowPattern.crossflow,
        "dynamic": False}))


all_elements = None
tube_in_element = None
tube_out_element = None
shell_in_elements = None

setup_block, arc_block = code_gen(n_passes, n_elements_per_pass, mixing=False)

exec(setup_block)

logger.info('Adding dynamic variables and constraints...')
for e in all_elements:
    e.setup()

logger.info('Applying time-discretization...')
m.discretizer = pe.TransformationFactory('dae.finite_difference')
m.discretizer.apply_to(m, nfe=100, wrt=m.fs.time, scheme="BACKWARD")

shell_inlet_temperature = 288.15
shell_flow = 44004.14222
tube_inlet_temperature = 384.35
tube_inlet_pressure = 7751362.5
tube_outlet_pressure = 7599375
tube_flow = 13896.84163

shell_area = 690073.9153
tube_area = 19542.2771

tube_mass = 1160 * 322

logger.info('Applying constraints...')
# Fix the heat transfer parameters in every element
for e in all_elements:
    e.crossflow_factor.fix(1)
    e.area.fix(1)
    e.tube_length = 195 / n_elements
    e.internal_surface_area = tube_area / n_elements
    e.external_surface_area = shell_area / n_elements
    e.heat_capacity = tube_mass * 466 / n_elements

tube_inlet_enthalpy = swco2.htpx(T=tube_inlet_temperature * pe.units.K, P=tube_inlet_pressure * pe.units.Pa)


def fix_tube_inlet(e):
    e.tube_inlet.flow_mol[:].fix(tube_flow)
    e.tube_inlet.pressure[:].fix(tube_inlet_pressure)
    e.tube_inlet.enth_mol[:].fix(tube_inlet_enthalpy)


def unfix_tube_inlet(e):
    e.tube_inlet.flow_mol[:].unfix()
    e.tube_inlet.pressure[:].unfix()
    e.tube_inlet.enth_mol[:].unfix()


def fix_shell_inlet(e):
    e.shell_inlet.flow_mol_comp[:, "H2O"].fix(0.01027 * shell_flow / n_elements_per_pass)
    e.shell_inlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * shell_flow / n_elements_per_pass)
    e.shell_inlet.flow_mol_comp[:, "N2"].fix(0.780066026 * shell_flow / n_elements_per_pass)
    e.shell_inlet.flow_mol_comp[:, "O2"].fix(0.209252382 * shell_flow / n_elements_per_pass)
    e.shell_inlet.flow_mol_comp[:, "NO"].fix(0)
    e.shell_inlet.flow_mol_comp[:, "SO2"].fix(0)
    e.shell_inlet.temperature[:].fix(shell_inlet_temperature)
    e.shell_inlet.pressure[:].fix(101325)


def unfix_shell_inlet(e):
    e.shell_inlet.flow_mol_comp[:, "H2O"].unfix()
    e.shell_inlet.flow_mol_comp[:, "CO2"].unfix()
    e.shell_inlet.flow_mol_comp[:, "N2"].unfix()
    e.shell_inlet.flow_mol_comp[:, "O2"].unfix()
    e.shell_inlet.flow_mol_comp[:, "NO"].unfix()
    e.shell_inlet.flow_mol_comp[:, "SO2"].unfix()
    e.shell_inlet.temperature[:].unfix()
    e.shell_inlet.pressure[:].unfix()


logger.info('Initializing elements...')

# Fix, initialize, then unfix
for e in all_elements:
    fix_shell_inlet(e)
    fix_tube_inlet(e)
    e.initialize()
    unfix_shell_inlet(e)
    unfix_tube_inlet(e)

# Reset inlet constraints
fix_tube_inlet(tube_in_element)
for e in shell_in_elements:
    fix_shell_inlet(e)


logger.info('Adding Arcs for 2D flow network...')

exec(arc_block)

# Apply Arc constraints
pe.TransformationFactory("network.expand_arcs").apply_to(m.fs)


solver = pe.SolverFactory('ipopt')
solver.options = {
            "tol": 1e-5,
            "linear_solver": "ma27",
            "max_iter": 500,
        }
logger.info('Solving model with steady state conditions...')
solver.solve(m, tee=True)

# Adding temperature disturbances
for e in shell_in_elements:
    for t in m.fs.time:
        if t >= 300 and t < 600:
            e.shell_inlet.temperature[t].fix(288.15 - 10)
        elif t >= 600 and t < 900:
            e.shell_inlet.temperature[t].fix(288.15)
        elif t >=900 and t < 1200:
            e.shell_inlet.temperature[t].fix(288.15 + 10)
        elif t >= 1200:
            e.shell_inlet.temperature[t].fix(288.15)


logger.info('Solving model with temperature step change...')
solver.solve(m, tee=True)


logger.info('Displaying results...')

hd = 0
for idx, e in enumerate(all_elements):
    print(f'')
    print(f'-------ELEMENT #{idx}-(t=0)------')
    hd += print_results_0d(e, t=0)

print('')
print(f'Total heat duty (t=0): {hd}')

hd = 0
for idx, e in enumerate(all_elements):
    print(f'')
    print(f'-------ELEMENT #{idx}-(t=600)------')
    hd += print_results_0d(e, t=600)

print('')
print(f'Total heat duty (t=600): {hd}')

res = []
for e in all_elements:
    t_tube_in = pe.value(e.tube.properties_in[0].temperature)
    t_tube_out = pe.value(e.tube.properties_out[0].temperature)
    t_shell_in = pe.value(e.shell.properties_in[0].temperature)
    t_shell_out = pe.value(e.shell.properties_out[0].temperature)
    res.append([t_tube_in, t_tube_out, t_shell_in, t_shell_out])

x = np.array(res)


fig, ax = plt.subplots(1, figsize=(12, 6))

for i in range(n_passes):
    ax.axvline(i * n_elements_per_pass, linestyle='dashed', c='grey', linewidth=1)

ax.plot(x.T[0], label='T tube in')
ax.plot(x.T[1], label='T tube out')
ax.plot(x.T[2], '.', label='T shell in')
ax.plot(x.T[3], '.', label='T shell out')

for i, t_in in enumerate(x.T[2]):
    y = (t_in, x.T[3][i])
    X = (i, i)
    ax.plot(X, y, '--', c='k', linewidth=0.5)

ax.legend()
ax.set_xlabel('Finite Element #')
ax.set_ylabel('Temperature (K)')
ax.set_title('2D Model Temperature Profile (t=0)')


res = []
for e in all_elements:
    t_tube_in = pe.value(e.tube.properties_in[600].temperature)
    t_tube_out = pe.value(e.tube.properties_out[600].temperature)
    t_shell_in = pe.value(e.shell.properties_in[600].temperature)
    t_shell_out = pe.value(e.shell.properties_out[600].temperature)
    res.append([t_tube_in, t_tube_out, t_shell_in, t_shell_out])

x = np.array(res)


fig, ax = plt.subplots(1, figsize=(12, 6))

for i in range(n_passes):
    ax.axvline(i * n_elements_per_pass, linestyle='dashed', c='grey', linewidth=1)

ax.plot(x.T[0], label='T tube in')
ax.plot(x.T[1], label='T tube out')
ax.plot(x.T[2], '.', label='T shell in')
ax.plot(x.T[3], '.', label='T shell out')

for i, t_in in enumerate(x.T[2]):
    y = (t_in, x.T[3][i])
    X = (i, i)
    ax.plot(X, y, '--', c='k', linewidth=0.5)

ax.legend()
ax.set_xlabel('Finite Element #')
ax.set_ylabel('Temperature (K)')
ax.set_title('2D Model Temperature Profile (t=600)')


t_out = pe.value(tube_out_element.tube.properties_out[:].temperature)
time = [t for t in m.fs.time]

fig, ax = plt.subplots(1)
ax.plot(time, t_out)
ax.set_xlabel('Time (s)')
ax.set_ylabel('CO2 Exit Temperature (K)')

plt.show()
