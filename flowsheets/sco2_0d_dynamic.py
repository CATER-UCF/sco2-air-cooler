"""
This a simple flowsheet containing a 0D, air-cooled heat exchanger. The unit
model is setup differently in four examples. Our ultimate goal is to model
countercurrent, crossflow geometry in two dimensions as a network of 0D models.
For now though, I just want to make sure I understand how transient modelling
works in IDAES.


Example 1: steady state
- This is a baseline case to get us started and compare.

Example 2: time discretized
- Here, the dynamic flag is set to False but the model is time discretized. My next step
  will be to add wall temperature and time-derivative terms for a lumped capacitance model.

Example 3: time discretized dynamic
- Here, the dynamic flag is set to True. The model solves but the results don't
  make sense since exit mass flows are free to change.

Example 4: time discretized dynamic fixed mass
- Same as Example 3 but the mass flows are all fixed. Now there are too few degrees of freedom.
"""
import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pe
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import HeatExchanger, HeatExchangerFlowPattern
from idaes.generic_models.unit_models.heat_exchanger import delta_temperature_amtd_callback
from idaes.generic_models.properties import swco2
from idaes.power_generation.properties import FlueGasParameterBlock
from util import print_results_0d


def set_boundary_conditions(m):

    shell_inlet_temperature = 288.15
    shell_flow = 44004.14222
    tube_inlet_temperature = 384.35
    tube_inlet_pressure = 7751362.5
    tube_outlet_pressure = 7599375
    tube_flow = 13896.84163

    shell_area = 690073.9153
    shell_HTC = 30

    tube_area = 19542.2771
    tube_HTC = 1000

    UA = 1 / (1 / (shell_area * shell_HTC) + 1 / (tube_area * tube_HTC))

    m.fs.HE.crossflow_factor.fix(0.8)
    m.fs.HE.area.fix(1)
    m.fs.HE.overall_heat_transfer_coefficient[:].fix(UA)

    tube_inlet_enthalpy = swco2.htpx(T=tube_inlet_temperature * pe.units.K, P=tube_inlet_pressure * pe.units.Pa)

    m.fs.HE.tube_inlet.flow_mol[:].fix(tube_flow)
    m.fs.HE.tube_inlet.pressure[:].fix(tube_inlet_pressure)
    m.fs.HE.tube_inlet.enth_mol[:].fix(tube_inlet_enthalpy)
    m.fs.HE.tube_outlet.pressure[:].fix(tube_outlet_pressure)

    m.fs.HE.shell_inlet.flow_mol_comp[:, "H2O"].fix(0.01027 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "N2"].fix(0.780066026 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "O2"].fix(0.209252382 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "NO"].fix(0)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "SO2"].fix(0)
    m.fs.HE.shell_inlet.temperature[:].fix(shell_inlet_temperature)
    m.fs.HE.shell_inlet.pressure[:].fix(101325)
    m.fs.HE.shell_outlet.pressure[:].fix(101325 * 0.95)
    return m


def fix_outlet_flows(m):

    shell_flow = 44004.14222
    tube_flow = 13896.84163

    m.fs.HE.tube_outlet.flow_mol[:].fix(tube_flow)

    m.fs.HE.shell_outlet.flow_mol_comp[:, "H2O"].fix(0.01027 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "N2"].fix(0.780066026 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "O2"].fix(0.209252382 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "NO"].fix(0)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "SO2"].fix(0)

    return m


def make_model(time_discretize=False, dyn=False):

    m = pe.ConcreteModel()
    if time_discretize:
        m.fs = FlowsheetBlock(default={"dynamic": True, "time_set": [0, 120], "time_units": pe.units.s})
    else:
        m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.HE = HeatExchanger(default={
        "delta_temperature_callback": delta_temperature_amtd_callback,
        "cold_side_name": "shell",
        "hot_side_name": "tube",
        "shell": {"property_package": m.fs.prop_fluegas,
                  "has_pressure_change": True},
        "tube": {"property_package": m.fs.prop_sco2,
                 "has_pressure_change": True},
        "flow_pattern": HeatExchangerFlowPattern.crossflow,
        "dynamic": dyn})

    if time_discretize:
        m.discretizer = pe.TransformationFactory('dae.finite_difference')
        m.discretizer.apply_to(m, nfe=100, wrt=m.fs.time, scheme="BACKWARD")

    return set_boundary_conditions(m)


"""
Example 1: Steady-state
"""
m = make_model()
m.fs.model_check()
print('Initializing steady-state model...')
m.fs.HE.initialize()

solver = pe.SolverFactory('ipopt')
solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 500,
        }
solver.solve(m, tee=True)

print('')
print('Steady-state results:')
print('+++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE)


"""
Example 2: Time-discretized model with temperature step change
"""
m = make_model(time_discretize=True)
m.fs.model_check()
print('Initializing time-discretized model...')
m.fs.HE.initialize()

solver = pe.SolverFactory('ipopt')
solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 500,
        }
solver.solve(m, tee=True)

for t in m.fs.time:
    if t > 60:
        m.fs.HE.shell_inlet.temperature[t].fix(288.15 + 10)

solver.solve(m, tee=True)

print('')
print('Time-discretized results (at time=0):')
print('+++++++++++++++++++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE)

print('')
print('Time-discretized results (at time=120):')
print('+++++++++++++++++++++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE, t=120)


"""
Example 3: Time-discretized dynamic model (no step change) 
"""
m = make_model(time_discretize=True, dyn=True)
m.fs.model_check()
print('Initializing dynamic model...')
m.fs.HE.initialize()

solver = pe.SolverFactory('ipopt')
solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 500,
        }
solver.solve(m, tee=True)

print('')
print('Dynamic model results (at time=0):')
print('++++++++++++++++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE)

print('')
print('Dynamic model results (at time=120):')
print('++++++++++++++++++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE, t=120)

t_tube_in = pe.value(m.fs.HE.tube.properties_in[:].temperature)
t_tube_out = pe.value(m.fs.HE.tube.properties_out[:].temperature)
t_shell_in = pe.value(m.fs.HE.shell.properties_in[:].temperature)
t_shell_out = pe.value(m.fs.HE.shell.properties_out[:].temperature)
hd = np.array(pe.value(m.fs.HE.heat_duty[:])) * -1e-6

g_tube_in = pe.value(m.fs.HE.tube.properties_in[:].flow_mol)
g_tube_out = pe.value(m.fs.HE.tube.properties_out[:].flow_mol)

plt.plot(t_tube_in, label='t tube in')
plt.plot(t_tube_out, label='t tube out')
plt.plot(t_shell_in, label='t shell in')
plt.plot(t_shell_out, label='t shell out')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Dynamic Model - Temperatures')
plt.legend()

fig = plt.subplots(1)
plt.plot(g_tube_in, label='Tube in')
plt.plot(g_tube_out, label='Tube out')
plt.xlabel('Time (s)')
plt.ylabel('Mass flow (mol/s)')
plt.title('Dynamic Model - sCO2 Mass Flows')
plt.legend()


h2o = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["H2O"])
co2 = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["CO2"])
n2 = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["N2"])
o2 = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["O2"])

fig = plt.subplots(1)
plt.plot(h2o, label='H2O out')
plt.plot(co2, label='CO2 out')
plt.plot(n2, label='N2 out')
plt.plot(o2, label='O2 out')
plt.xlabel('Time (s)')
plt.ylabel('Mass flow (mol/s)')
plt.title('Dynamic Model - Air Component Mass Flows')
plt.legend()

plt.show()


"""
Example 4: Dynamic model with fixed mass flows 
"""
m = make_model(time_discretize=True, dyn=True)
m = fix_outlet_flows(m)
m.fs.model_check()
print('Initializing dynamic model...')
m.fs.HE.initialize()

solver = pe.SolverFactory('ipopt')
solver.options = {
            "tol": 1e-7,
            "linear_solver": "ma27",
            "max_iter": 500,
        }
solver.solve(m, tee=True)
