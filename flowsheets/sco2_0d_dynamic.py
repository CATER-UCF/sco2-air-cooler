"""
Work in progress...
"""
import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pe
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import HeatExchanger, HeatExchangerFlowPattern
from idaes.generic_models.unit_models.heat_exchanger import delta_temperature_amtd_callback
from idaes.generic_models.properties import swco2
from idaes.power_generation.properties import FlueGasParameterBlock
from util import print_results_0d

#from models import HeatExchangerDynamic


def set_boundary_conditions(m):

    shell_inlet_temperature = 288.15
    shell_flow = 44004.14222
    tube_inlet_temperature = 384.35
    tube_inlet_pressure = 7751362.5
    tube_outlet_pressure = 7599375
    tube_flow = 13896.84163

    shell_area = 690073.9153
    shell_HTC = 100

    tube_area = 19542.2771
    tube_HTC = 2000

    UA = 1 / (1 / (shell_area * shell_HTC) + 1 / (tube_area * tube_HTC))

    m.fs.HE.crossflow_factor.fix(0.9)
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

    """
    m.fs.HE.tube_outlet.flow_mol[:].fix(tube_flow)
    m.fs.HE.tube_outlet.flow_mol[:].unfix()

    m.fs.HE.shell_outlet.flow_mol_comp[:, "H2O"].fix(0.01027 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "N2"].fix(0.780066026 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "O2"].fix(0.209252382 * shell_flow)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "NO"].fix(0)
    m.fs.HE.shell_outlet.flow_mol_comp[:, "SO2"].fix(0)

    m.fs.HE.shell_outlet.flow_mol_comp[:, "H2O"].unfix()
    m.fs.HE.shell_outlet.flow_mol_comp[:, "CO2"].unfix()
    m.fs.HE.shell_outlet.flow_mol_comp[:, "N2"].unfix()
    m.fs.HE.shell_outlet.flow_mol_comp[:, "O2"].unfix()
    m.fs.HE.shell_outlet.flow_mol_comp[:, "NO"].unfix()
    m.fs.HE.shell_outlet.flow_mol_comp[:, "SO2"].unfix()
    """

    return m


def make_model(dyn=False):

    m = pe.ConcreteModel()
    if dyn:
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

    if dyn:
        m.discretizer = pe.TransformationFactory('dae.finite_difference')
        m.discretizer.apply_to(m, nfe=10, wrt=m.fs.time, scheme="BACKWARD")

    return set_boundary_conditions(m)


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


m = make_model(dyn=True)
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
print('Dynamic results (at time=0):')
print('++++++++++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE)


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
plt.legend()

fig = plt.subplots(1)
plt.plot(g_tube_in, label='g tube in')
plt.plot(g_tube_out, label='g tube out')
plt.legend()


h2o = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["H2O"])
co2 = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["CO2"])
n2 = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["N2"])
o2 = pe.value(m.fs.HE.shell.properties_out[:].flow_mol_comp["O2"])

fig = plt.subplots(1)
plt.plot(h2o, label='h2o shell out')
plt.plot(co2, label='co2 shell out')
plt.plot(n2, label='n2 shell out')
plt.plot(o2, label='o2 shell out')
plt.legend()

plt.show()
