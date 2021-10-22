import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pe
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import HeatExchanger, HeatExchangerFlowPattern
from idaes.generic_models.unit_models.heat_exchanger import delta_temperature_lmtd_callback
from idaes.generic_models.properties import swco2
from idaes.power_generation.properties import FlueGasParameterBlock
from util import print_results_0d, write_csv
from models import HeatExchangerElement


def set_boundary_conditions(m):

    shell_inlet_temperature = 288.15
    shell_flow = 44004.14222
    tube_inlet_temperature = 384.35
    tube_inlet_pressure = 7653000
    tube_flow = 13896.84163

    shell_area = 690073.9153
    tube_area = 19542.2771
    tube_mass = 1160 * 322

    m.fs.HE.tube_length = 195
    m.fs.HE.internal_surface_area = tube_area
    m.fs.HE.external_surface_area = shell_area

    m.fs.HE.crossflow_factor.fix(0.8)
    m.fs.HE.area.fix(1)
    m.fs.HE.heat_capacity = tube_mass * 466

    tube_inlet_enthalpy = swco2.htpx(T=tube_inlet_temperature * pe.units.K, P=tube_inlet_pressure * pe.units.Pa)

    m.fs.HE.tube_inlet.flow_mol[:].fix(tube_flow)
    m.fs.HE.tube_inlet.pressure[:].fix(tube_inlet_pressure)
    m.fs.HE.tube_inlet.enth_mol[:].fix(tube_inlet_enthalpy)

    m.fs.HE.shell_inlet.flow_mol_comp[:, "H2O"].fix(0.01027 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "N2"].fix(0.780066026 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "O2"].fix(0.209252382 * shell_flow)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "NO"].fix(0)
    m.fs.HE.shell_inlet.flow_mol_comp[:, "SO2"].fix(0)
    m.fs.HE.shell_inlet.temperature[:].fix(shell_inlet_temperature)
    m.fs.HE.shell_inlet.pressure[:].fix(101325)

    m.fs.HE.tube_outlet[:].enth_mol.setub(tube_inlet_enthalpy)
    m.fs.HE.shell_outlet.temperature[:].setlb(shell_inlet_temperature)

    return m


def make_model():

    m = pe.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   "time_set": [0, 300, 600, 900, 1200, 1500],
                                   "time_units": pe.units.s})

    m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.HE = HeatExchangerElement(default={
        "delta_temperature_callback": delta_temperature_lmtd_callback,
        "cold_side_name": "shell",
        "hot_side_name": "tube",
        "shell": {"property_package": m.fs.prop_fluegas,
                  "has_pressure_change": False},
        "tube": {"property_package": m.fs.prop_sco2,
                 "has_pressure_change": False},
        "flow_pattern": HeatExchangerFlowPattern.crossflow,
        "dynamic": False})

    m.fs.HE.setup()
    m.fs.HE.activate_dynamic_heat_eq()
    m.discretizer = pe.TransformationFactory('dae.finite_difference')
    m.discretizer.apply_to(m, nfe=200, wrt=m.fs.time, scheme="BACKWARD")

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
            "tol": 1e-6,
            "linear_solver": "ma27",
            "max_iter": 500,
        }
solver.solve(m, tee=True)

print('')
print('Steady-state results:')
print('+++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE)

t_wall = pe.value(m.fs.HE.wall_temperature[0])
print(f'T wall: {t_wall}')

for t in m.fs.time:
    if t >= 300 and t < 600:
        m.fs.HE.shell_inlet.temperature[t].fix(288.15 - 10)
    elif t >= 600 and t < 900:
        m.fs.HE.shell_inlet.temperature[t].fix(288.15)
    elif t >=900 and t < 1200:
        m.fs.HE.shell_inlet.temperature[t].fix(288.15 + 10)
    elif t >= 1200:
        m.fs.HE.shell_inlet.temperature[t].fix(288.15)

solver.solve(m, tee=True)

print('')
print('Time-discretized results (at time=0):')
print('+++++++++++++++++++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE)

print('')
print('Time-discretized results (at time=600):')
print('+++++++++++++++++++++++++++++++++++++++')
print('')
print_results_0d(m.fs.HE, t=600)

write_csv(f'./data/time_series_0d.csv', [m.fs.HE])
