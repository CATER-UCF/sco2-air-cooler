"""
Steady-state simulation of a 2D countercurrent crossflow heat exchanger.
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
from util import print_results_0d

m = pe.ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
m.fs.prop_fluegas = FlueGasParameterBlock()

n_passes = 8
n_elements_per_pass = 7
n_elements = n_passes * n_elements_per_pass
m.fs.es = []

for _ in range(n_elements):
    m.fs.es.append(HeatExchanger(default={
        "delta_temperature_callback": delta_temperature_lmtd_callback,
        "cold_side_name": "shell",
        "hot_side_name": "tube",
        "shell": {"property_package": m.fs.prop_fluegas,
                  "has_pressure_change": True},
        "tube": {"property_package": m.fs.prop_sco2,
                "has_pressure_change": True},
        "flow_pattern": HeatExchangerFlowPattern.crossflow}))


# Paste codgen setup code block here:
m.fs.e0 = m.fs.es[0]
m.fs.e1 = m.fs.es[1]
m.fs.e2 = m.fs.es[2]
m.fs.e3 = m.fs.es[3]
m.fs.e4 = m.fs.es[4]
m.fs.e5 = m.fs.es[5]
m.fs.e6 = m.fs.es[6]
m.fs.e7 = m.fs.es[7]
m.fs.e8 = m.fs.es[8]
m.fs.e9 = m.fs.es[9]
m.fs.e10 = m.fs.es[10]
m.fs.e11 = m.fs.es[11]
m.fs.e12 = m.fs.es[12]
m.fs.e13 = m.fs.es[13]
m.fs.e14 = m.fs.es[14]
m.fs.e15 = m.fs.es[15]
m.fs.e16 = m.fs.es[16]
m.fs.e17 = m.fs.es[17]
m.fs.e18 = m.fs.es[18]
m.fs.e19 = m.fs.es[19]
m.fs.e20 = m.fs.es[20]
m.fs.e21 = m.fs.es[21]
m.fs.e22 = m.fs.es[22]
m.fs.e23 = m.fs.es[23]
m.fs.e24 = m.fs.es[24]
m.fs.e25 = m.fs.es[25]
m.fs.e26 = m.fs.es[26]
m.fs.e27 = m.fs.es[27]
m.fs.e28 = m.fs.es[28]
m.fs.e29 = m.fs.es[29]
m.fs.e30 = m.fs.es[30]
m.fs.e31 = m.fs.es[31]
m.fs.e32 = m.fs.es[32]
m.fs.e33 = m.fs.es[33]
m.fs.e34 = m.fs.es[34]
m.fs.e35 = m.fs.es[35]
m.fs.e36 = m.fs.es[36]
m.fs.e37 = m.fs.es[37]
m.fs.e38 = m.fs.es[38]
m.fs.e39 = m.fs.es[39]
m.fs.e40 = m.fs.es[40]
m.fs.e41 = m.fs.es[41]
m.fs.e42 = m.fs.es[42]
m.fs.e43 = m.fs.es[43]
m.fs.e44 = m.fs.es[44]
m.fs.e45 = m.fs.es[45]
m.fs.e46 = m.fs.es[46]
m.fs.e47 = m.fs.es[47]
m.fs.e48 = m.fs.es[48]
m.fs.e49 = m.fs.es[49]
m.fs.e50 = m.fs.es[50]
m.fs.e51 = m.fs.es[51]
m.fs.e52 = m.fs.es[52]
m.fs.e53 = m.fs.es[53]
m.fs.e54 = m.fs.es[54]
m.fs.e55 = m.fs.es[55]

all_elements = [m.fs.e0, m.fs.e1, m.fs.e2, m.fs.e3, m.fs.e4, m.fs.e5, m.fs.e6,
    m.fs.e7, m.fs.e8, m.fs.e9, m.fs.e10, m.fs.e11, m.fs.e12, m.fs.e13,
    m.fs.e14, m.fs.e15, m.fs.e16, m.fs.e17, m.fs.e18, m.fs.e19, m.fs.e20,
    m.fs.e21, m.fs.e22, m.fs.e23, m.fs.e24, m.fs.e25, m.fs.e26, m.fs.e27,
    m.fs.e28, m.fs.e29, m.fs.e30, m.fs.e31, m.fs.e32, m.fs.e33, m.fs.e34,
    m.fs.e35, m.fs.e36, m.fs.e37, m.fs.e38, m.fs.e39, m.fs.e40, m.fs.e41,
    m.fs.e42, m.fs.e43, m.fs.e44, m.fs.e45, m.fs.e46, m.fs.e47, m.fs.e48,
    m.fs.e49, m.fs.e50, m.fs.e51, m.fs.e52, m.fs.e53, m.fs.e54, m.fs.e55]

shell_in_elements = [m.fs.e49, m.fs.e50, m.fs.e51, m.fs.e52, m.fs.e53,
    m.fs.e54, m.fs.e55]

tube_in_element = m.fs.e0
tube_out_element = m.fs.e55

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

# Fix the HTCs in every element
for e in all_elements:
    e.crossflow_factor.fix(0.5)
    e.area.fix(1)
    e.overall_heat_transfer_coefficient[:].fix(UA / n_elements)

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
    e.shell.deltaP[:].fix(0)


def unfix_shell_inlet(e):
    e.shell_inlet.flow_mol_comp[:, "H2O"].unfix()
    e.shell_inlet.flow_mol_comp[:, "CO2"].unfix()
    e.shell_inlet.flow_mol_comp[:, "N2"].unfix()
    e.shell_inlet.flow_mol_comp[:, "O2"].unfix()
    e.shell_inlet.flow_mol_comp[:, "NO"].unfix()
    e.shell_inlet.flow_mol_comp[:, "SO2"].unfix()
    e.shell_inlet.temperature[:].unfix()
    e.shell_inlet.pressure[:].unfix()


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

# Tube-side pressure loss
tube_dp = (tube_outlet_pressure - tube_inlet_pressure) / n_elements
tube_p = tube_inlet_pressure
for idx, e in enumerate(all_elements):
    tube_p += tube_dp
    e.tube_outlet.pressure[:].fix(tube_p)

# Paste codgen arc code block here:

# Add tube-side Arcs
m.fs.e_tube_Arc0 = Arc(source=m.fs.e0.outlet_1, destination=m.fs.e1.inlet_1)
m.fs.e_tube_Arc1 = Arc(source=m.fs.e1.outlet_1, destination=m.fs.e2.inlet_1)
m.fs.e_tube_Arc2 = Arc(source=m.fs.e2.outlet_1, destination=m.fs.e3.inlet_1)
m.fs.e_tube_Arc3 = Arc(source=m.fs.e3.outlet_1, destination=m.fs.e4.inlet_1)
m.fs.e_tube_Arc4 = Arc(source=m.fs.e4.outlet_1, destination=m.fs.e5.inlet_1)
m.fs.e_tube_Arc5 = Arc(source=m.fs.e5.outlet_1, destination=m.fs.e6.inlet_1)
m.fs.e_tube_Arc6 = Arc(source=m.fs.e6.outlet_1, destination=m.fs.e7.inlet_1)
m.fs.e_tube_Arc7 = Arc(source=m.fs.e7.outlet_1, destination=m.fs.e8.inlet_1)
m.fs.e_tube_Arc8 = Arc(source=m.fs.e8.outlet_1, destination=m.fs.e9.inlet_1)
m.fs.e_tube_Arc9 = Arc(source=m.fs.e9.outlet_1, destination=m.fs.e10.inlet_1)
m.fs.e_tube_Arc10 = Arc(source=m.fs.e10.outlet_1, destination=m.fs.e11.inlet_1)
m.fs.e_tube_Arc11 = Arc(source=m.fs.e11.outlet_1, destination=m.fs.e12.inlet_1)
m.fs.e_tube_Arc12 = Arc(source=m.fs.e12.outlet_1, destination=m.fs.e13.inlet_1)
m.fs.e_tube_Arc13 = Arc(source=m.fs.e13.outlet_1, destination=m.fs.e14.inlet_1)
m.fs.e_tube_Arc14 = Arc(source=m.fs.e14.outlet_1, destination=m.fs.e15.inlet_1)
m.fs.e_tube_Arc15 = Arc(source=m.fs.e15.outlet_1, destination=m.fs.e16.inlet_1)
m.fs.e_tube_Arc16 = Arc(source=m.fs.e16.outlet_1, destination=m.fs.e17.inlet_1)
m.fs.e_tube_Arc17 = Arc(source=m.fs.e17.outlet_1, destination=m.fs.e18.inlet_1)
m.fs.e_tube_Arc18 = Arc(source=m.fs.e18.outlet_1, destination=m.fs.e19.inlet_1)
m.fs.e_tube_Arc19 = Arc(source=m.fs.e19.outlet_1, destination=m.fs.e20.inlet_1)
m.fs.e_tube_Arc20 = Arc(source=m.fs.e20.outlet_1, destination=m.fs.e21.inlet_1)
m.fs.e_tube_Arc21 = Arc(source=m.fs.e21.outlet_1, destination=m.fs.e22.inlet_1)
m.fs.e_tube_Arc22 = Arc(source=m.fs.e22.outlet_1, destination=m.fs.e23.inlet_1)
m.fs.e_tube_Arc23 = Arc(source=m.fs.e23.outlet_1, destination=m.fs.e24.inlet_1)
m.fs.e_tube_Arc24 = Arc(source=m.fs.e24.outlet_1, destination=m.fs.e25.inlet_1)
m.fs.e_tube_Arc25 = Arc(source=m.fs.e25.outlet_1, destination=m.fs.e26.inlet_1)
m.fs.e_tube_Arc26 = Arc(source=m.fs.e26.outlet_1, destination=m.fs.e27.inlet_1)
m.fs.e_tube_Arc27 = Arc(source=m.fs.e27.outlet_1, destination=m.fs.e28.inlet_1)
m.fs.e_tube_Arc28 = Arc(source=m.fs.e28.outlet_1, destination=m.fs.e29.inlet_1)
m.fs.e_tube_Arc29 = Arc(source=m.fs.e29.outlet_1, destination=m.fs.e30.inlet_1)
m.fs.e_tube_Arc30 = Arc(source=m.fs.e30.outlet_1, destination=m.fs.e31.inlet_1)
m.fs.e_tube_Arc31 = Arc(source=m.fs.e31.outlet_1, destination=m.fs.e32.inlet_1)
m.fs.e_tube_Arc32 = Arc(source=m.fs.e32.outlet_1, destination=m.fs.e33.inlet_1)
m.fs.e_tube_Arc33 = Arc(source=m.fs.e33.outlet_1, destination=m.fs.e34.inlet_1)
m.fs.e_tube_Arc34 = Arc(source=m.fs.e34.outlet_1, destination=m.fs.e35.inlet_1)
m.fs.e_tube_Arc35 = Arc(source=m.fs.e35.outlet_1, destination=m.fs.e36.inlet_1)
m.fs.e_tube_Arc36 = Arc(source=m.fs.e36.outlet_1, destination=m.fs.e37.inlet_1)
m.fs.e_tube_Arc37 = Arc(source=m.fs.e37.outlet_1, destination=m.fs.e38.inlet_1)
m.fs.e_tube_Arc38 = Arc(source=m.fs.e38.outlet_1, destination=m.fs.e39.inlet_1)
m.fs.e_tube_Arc39 = Arc(source=m.fs.e39.outlet_1, destination=m.fs.e40.inlet_1)
m.fs.e_tube_Arc40 = Arc(source=m.fs.e40.outlet_1, destination=m.fs.e41.inlet_1)
m.fs.e_tube_Arc41 = Arc(source=m.fs.e41.outlet_1, destination=m.fs.e42.inlet_1)
m.fs.e_tube_Arc42 = Arc(source=m.fs.e42.outlet_1, destination=m.fs.e43.inlet_1)
m.fs.e_tube_Arc43 = Arc(source=m.fs.e43.outlet_1, destination=m.fs.e44.inlet_1)
m.fs.e_tube_Arc44 = Arc(source=m.fs.e44.outlet_1, destination=m.fs.e45.inlet_1)
m.fs.e_tube_Arc45 = Arc(source=m.fs.e45.outlet_1, destination=m.fs.e46.inlet_1)
m.fs.e_tube_Arc46 = Arc(source=m.fs.e46.outlet_1, destination=m.fs.e47.inlet_1)
m.fs.e_tube_Arc47 = Arc(source=m.fs.e47.outlet_1, destination=m.fs.e48.inlet_1)
m.fs.e_tube_Arc48 = Arc(source=m.fs.e48.outlet_1, destination=m.fs.e49.inlet_1)
m.fs.e_tube_Arc49 = Arc(source=m.fs.e49.outlet_1, destination=m.fs.e50.inlet_1)
m.fs.e_tube_Arc50 = Arc(source=m.fs.e50.outlet_1, destination=m.fs.e51.inlet_1)
m.fs.e_tube_Arc51 = Arc(source=m.fs.e51.outlet_1, destination=m.fs.e52.inlet_1)
m.fs.e_tube_Arc52 = Arc(source=m.fs.e52.outlet_1, destination=m.fs.e53.inlet_1)
m.fs.e_tube_Arc53 = Arc(source=m.fs.e53.outlet_1, destination=m.fs.e54.inlet_1)
m.fs.e_tube_Arc54 = Arc(source=m.fs.e54.outlet_1, destination=m.fs.e55.inlet_1)

# Shell mixer & separator setup
m.fs.mix0 = Mixer(default={"property_package": m.fs.prop_fluegas, "num_inlets": 7})
m.fs.mix1 = Mixer(default={"property_package": m.fs.prop_fluegas, "num_inlets": 7})
m.fs.mix2 = Mixer(default={"property_package": m.fs.prop_fluegas, "num_inlets": 7})
m.fs.mix3 = Mixer(default={"property_package": m.fs.prop_fluegas, "num_inlets": 7})
m.fs.mix4 = Mixer(default={"property_package": m.fs.prop_fluegas, "num_inlets": 7})
m.fs.mix5 = Mixer(default={"property_package": m.fs.prop_fluegas, "num_inlets": 7})
m.fs.mix6 = Mixer(default={"property_package": m.fs.prop_fluegas, "num_inlets": 7})

m.fs.sep0 = Separator(default={"property_package": m.fs.prop_fluegas, "num_outlets": 7})
m.fs.sep1 = Separator(default={"property_package": m.fs.prop_fluegas, "num_outlets": 7})
m.fs.sep2 = Separator(default={"property_package": m.fs.prop_fluegas, "num_outlets": 7})
m.fs.sep3 = Separator(default={"property_package": m.fs.prop_fluegas, "num_outlets": 7})
m.fs.sep4 = Separator(default={"property_package": m.fs.prop_fluegas, "num_outlets": 7})
m.fs.sep5 = Separator(default={"property_package": m.fs.prop_fluegas, "num_outlets": 7})
m.fs.sep6 = Separator(default={"property_package": m.fs.prop_fluegas, "num_outlets": 7})

# Mixer Arcs for pass 1
m.fs.e_shell_pass1_m0 = Arc(source=m.fs.e49.outlet_2, destination=m.fs.mix0.inlet_1)
m.fs.e_shell_pass1_m1 = Arc(source=m.fs.e50.outlet_2, destination=m.fs.mix0.inlet_2)
m.fs.e_shell_pass1_m2 = Arc(source=m.fs.e51.outlet_2, destination=m.fs.mix0.inlet_3)
m.fs.e_shell_pass1_m3 = Arc(source=m.fs.e52.outlet_2, destination=m.fs.mix0.inlet_4)
m.fs.e_shell_pass1_m4 = Arc(source=m.fs.e53.outlet_2, destination=m.fs.mix0.inlet_5)
m.fs.e_shell_pass1_m5 = Arc(source=m.fs.e54.outlet_2, destination=m.fs.mix0.inlet_6)
m.fs.e_shell_pass1_m6 = Arc(source=m.fs.e55.outlet_2, destination=m.fs.mix0.inlet_7)

# Mixer / Separator Arc for pass 1
m.fs.e_shell_pass1_ms_Arc = Arc(source=m.fs.mix0.outlet, destination=m.fs.sep0.inlet)

# Separator Arcs for pass 1
m.fs.e_shell_pass1_s0 = Arc(source=m.fs.sep0.outlet_1, destination=m.fs.e42.inlet_2)
m.fs.e_shell_pass1_s1 = Arc(source=m.fs.sep0.outlet_2, destination=m.fs.e43.inlet_2)
m.fs.e_shell_pass1_s2 = Arc(source=m.fs.sep0.outlet_3, destination=m.fs.e44.inlet_2)
m.fs.e_shell_pass1_s3 = Arc(source=m.fs.sep0.outlet_4, destination=m.fs.e45.inlet_2)
m.fs.e_shell_pass1_s4 = Arc(source=m.fs.sep0.outlet_5, destination=m.fs.e46.inlet_2)
m.fs.e_shell_pass1_s5 = Arc(source=m.fs.sep0.outlet_6, destination=m.fs.e47.inlet_2)
m.fs.e_shell_pass1_s6 = Arc(source=m.fs.sep0.outlet_7, destination=m.fs.e48.inlet_2)

# Mixer Arcs for pass 2
m.fs.e_shell_pass2_m0 = Arc(source=m.fs.e42.outlet_2, destination=m.fs.mix1.inlet_1)
m.fs.e_shell_pass2_m1 = Arc(source=m.fs.e43.outlet_2, destination=m.fs.mix1.inlet_2)
m.fs.e_shell_pass2_m2 = Arc(source=m.fs.e44.outlet_2, destination=m.fs.mix1.inlet_3)
m.fs.e_shell_pass2_m3 = Arc(source=m.fs.e45.outlet_2, destination=m.fs.mix1.inlet_4)
m.fs.e_shell_pass2_m4 = Arc(source=m.fs.e46.outlet_2, destination=m.fs.mix1.inlet_5)
m.fs.e_shell_pass2_m5 = Arc(source=m.fs.e47.outlet_2, destination=m.fs.mix1.inlet_6)
m.fs.e_shell_pass2_m6 = Arc(source=m.fs.e48.outlet_2, destination=m.fs.mix1.inlet_7)

# Mixer / Separator Arc for pass 2
m.fs.e_shell_pass2_ms_Arc = Arc(source=m.fs.mix1.outlet, destination=m.fs.sep1.inlet)

# Separator Arcs for pass 2
m.fs.e_shell_pass2_s0 = Arc(source=m.fs.sep1.outlet_1, destination=m.fs.e35.inlet_2)
m.fs.e_shell_pass2_s1 = Arc(source=m.fs.sep1.outlet_2, destination=m.fs.e36.inlet_2)
m.fs.e_shell_pass2_s2 = Arc(source=m.fs.sep1.outlet_3, destination=m.fs.e37.inlet_2)
m.fs.e_shell_pass2_s3 = Arc(source=m.fs.sep1.outlet_4, destination=m.fs.e38.inlet_2)
m.fs.e_shell_pass2_s4 = Arc(source=m.fs.sep1.outlet_5, destination=m.fs.e39.inlet_2)
m.fs.e_shell_pass2_s5 = Arc(source=m.fs.sep1.outlet_6, destination=m.fs.e40.inlet_2)
m.fs.e_shell_pass2_s6 = Arc(source=m.fs.sep1.outlet_7, destination=m.fs.e41.inlet_2)

# Mixer Arcs for pass 3
m.fs.e_shell_pass3_m0 = Arc(source=m.fs.e35.outlet_2, destination=m.fs.mix2.inlet_1)
m.fs.e_shell_pass3_m1 = Arc(source=m.fs.e36.outlet_2, destination=m.fs.mix2.inlet_2)
m.fs.e_shell_pass3_m2 = Arc(source=m.fs.e37.outlet_2, destination=m.fs.mix2.inlet_3)
m.fs.e_shell_pass3_m3 = Arc(source=m.fs.e38.outlet_2, destination=m.fs.mix2.inlet_4)
m.fs.e_shell_pass3_m4 = Arc(source=m.fs.e39.outlet_2, destination=m.fs.mix2.inlet_5)
m.fs.e_shell_pass3_m5 = Arc(source=m.fs.e40.outlet_2, destination=m.fs.mix2.inlet_6)
m.fs.e_shell_pass3_m6 = Arc(source=m.fs.e41.outlet_2, destination=m.fs.mix2.inlet_7)

# Mixer / Separator Arc for pass 3
m.fs.e_shell_pass3_ms_Arc = Arc(source=m.fs.mix2.outlet, destination=m.fs.sep2.inlet)

# Separator Arcs for pass 3
m.fs.e_shell_pass3_s0 = Arc(source=m.fs.sep2.outlet_1, destination=m.fs.e28.inlet_2)
m.fs.e_shell_pass3_s1 = Arc(source=m.fs.sep2.outlet_2, destination=m.fs.e29.inlet_2)
m.fs.e_shell_pass3_s2 = Arc(source=m.fs.sep2.outlet_3, destination=m.fs.e30.inlet_2)
m.fs.e_shell_pass3_s3 = Arc(source=m.fs.sep2.outlet_4, destination=m.fs.e31.inlet_2)
m.fs.e_shell_pass3_s4 = Arc(source=m.fs.sep2.outlet_5, destination=m.fs.e32.inlet_2)
m.fs.e_shell_pass3_s5 = Arc(source=m.fs.sep2.outlet_6, destination=m.fs.e33.inlet_2)
m.fs.e_shell_pass3_s6 = Arc(source=m.fs.sep2.outlet_7, destination=m.fs.e34.inlet_2)

# Mixer Arcs for pass 4
m.fs.e_shell_pass4_m0 = Arc(source=m.fs.e28.outlet_2, destination=m.fs.mix3.inlet_1)
m.fs.e_shell_pass4_m1 = Arc(source=m.fs.e29.outlet_2, destination=m.fs.mix3.inlet_2)
m.fs.e_shell_pass4_m2 = Arc(source=m.fs.e30.outlet_2, destination=m.fs.mix3.inlet_3)
m.fs.e_shell_pass4_m3 = Arc(source=m.fs.e31.outlet_2, destination=m.fs.mix3.inlet_4)
m.fs.e_shell_pass4_m4 = Arc(source=m.fs.e32.outlet_2, destination=m.fs.mix3.inlet_5)
m.fs.e_shell_pass4_m5 = Arc(source=m.fs.e33.outlet_2, destination=m.fs.mix3.inlet_6)
m.fs.e_shell_pass4_m6 = Arc(source=m.fs.e34.outlet_2, destination=m.fs.mix3.inlet_7)

# Mixer / Separator Arc for pass 4
m.fs.e_shell_pass4_ms_Arc = Arc(source=m.fs.mix3.outlet, destination=m.fs.sep3.inlet)

# Separator Arcs for pass 4
m.fs.e_shell_pass4_s0 = Arc(source=m.fs.sep3.outlet_1, destination=m.fs.e21.inlet_2)
m.fs.e_shell_pass4_s1 = Arc(source=m.fs.sep3.outlet_2, destination=m.fs.e22.inlet_2)
m.fs.e_shell_pass4_s2 = Arc(source=m.fs.sep3.outlet_3, destination=m.fs.e23.inlet_2)
m.fs.e_shell_pass4_s3 = Arc(source=m.fs.sep3.outlet_4, destination=m.fs.e24.inlet_2)
m.fs.e_shell_pass4_s4 = Arc(source=m.fs.sep3.outlet_5, destination=m.fs.e25.inlet_2)
m.fs.e_shell_pass4_s5 = Arc(source=m.fs.sep3.outlet_6, destination=m.fs.e26.inlet_2)
m.fs.e_shell_pass4_s6 = Arc(source=m.fs.sep3.outlet_7, destination=m.fs.e27.inlet_2)

# Mixer Arcs for pass 5
m.fs.e_shell_pass5_m0 = Arc(source=m.fs.e21.outlet_2, destination=m.fs.mix4.inlet_1)
m.fs.e_shell_pass5_m1 = Arc(source=m.fs.e22.outlet_2, destination=m.fs.mix4.inlet_2)
m.fs.e_shell_pass5_m2 = Arc(source=m.fs.e23.outlet_2, destination=m.fs.mix4.inlet_3)
m.fs.e_shell_pass5_m3 = Arc(source=m.fs.e24.outlet_2, destination=m.fs.mix4.inlet_4)
m.fs.e_shell_pass5_m4 = Arc(source=m.fs.e25.outlet_2, destination=m.fs.mix4.inlet_5)
m.fs.e_shell_pass5_m5 = Arc(source=m.fs.e26.outlet_2, destination=m.fs.mix4.inlet_6)
m.fs.e_shell_pass5_m6 = Arc(source=m.fs.e27.outlet_2, destination=m.fs.mix4.inlet_7)

# Mixer / Separator Arc for pass 5
m.fs.e_shell_pass5_ms_Arc = Arc(source=m.fs.mix4.outlet, destination=m.fs.sep4.inlet)

# Separator Arcs for pass 5
m.fs.e_shell_pass5_s0 = Arc(source=m.fs.sep4.outlet_1, destination=m.fs.e14.inlet_2)
m.fs.e_shell_pass5_s1 = Arc(source=m.fs.sep4.outlet_2, destination=m.fs.e15.inlet_2)
m.fs.e_shell_pass5_s2 = Arc(source=m.fs.sep4.outlet_3, destination=m.fs.e16.inlet_2)
m.fs.e_shell_pass5_s3 = Arc(source=m.fs.sep4.outlet_4, destination=m.fs.e17.inlet_2)
m.fs.e_shell_pass5_s4 = Arc(source=m.fs.sep4.outlet_5, destination=m.fs.e18.inlet_2)
m.fs.e_shell_pass5_s5 = Arc(source=m.fs.sep4.outlet_6, destination=m.fs.e19.inlet_2)
m.fs.e_shell_pass5_s6 = Arc(source=m.fs.sep4.outlet_7, destination=m.fs.e20.inlet_2)

# Mixer Arcs for pass 6
m.fs.e_shell_pass6_m0 = Arc(source=m.fs.e14.outlet_2, destination=m.fs.mix5.inlet_1)
m.fs.e_shell_pass6_m1 = Arc(source=m.fs.e15.outlet_2, destination=m.fs.mix5.inlet_2)
m.fs.e_shell_pass6_m2 = Arc(source=m.fs.e16.outlet_2, destination=m.fs.mix5.inlet_3)
m.fs.e_shell_pass6_m3 = Arc(source=m.fs.e17.outlet_2, destination=m.fs.mix5.inlet_4)
m.fs.e_shell_pass6_m4 = Arc(source=m.fs.e18.outlet_2, destination=m.fs.mix5.inlet_5)
m.fs.e_shell_pass6_m5 = Arc(source=m.fs.e19.outlet_2, destination=m.fs.mix5.inlet_6)
m.fs.e_shell_pass6_m6 = Arc(source=m.fs.e20.outlet_2, destination=m.fs.mix5.inlet_7)

# Mixer / Separator Arc for pass 6
m.fs.e_shell_pass6_ms_Arc = Arc(source=m.fs.mix5.outlet, destination=m.fs.sep5.inlet)

# Separator Arcs for pass 6
m.fs.e_shell_pass6_s0 = Arc(source=m.fs.sep5.outlet_1, destination=m.fs.e7.inlet_2)
m.fs.e_shell_pass6_s1 = Arc(source=m.fs.sep5.outlet_2, destination=m.fs.e8.inlet_2)
m.fs.e_shell_pass6_s2 = Arc(source=m.fs.sep5.outlet_3, destination=m.fs.e9.inlet_2)
m.fs.e_shell_pass6_s3 = Arc(source=m.fs.sep5.outlet_4, destination=m.fs.e10.inlet_2)
m.fs.e_shell_pass6_s4 = Arc(source=m.fs.sep5.outlet_5, destination=m.fs.e11.inlet_2)
m.fs.e_shell_pass6_s5 = Arc(source=m.fs.sep5.outlet_6, destination=m.fs.e12.inlet_2)
m.fs.e_shell_pass6_s6 = Arc(source=m.fs.sep5.outlet_7, destination=m.fs.e13.inlet_2)

# Mixer Arcs for pass 7
m.fs.e_shell_pass7_m0 = Arc(source=m.fs.e7.outlet_2, destination=m.fs.mix6.inlet_1)
m.fs.e_shell_pass7_m1 = Arc(source=m.fs.e8.outlet_2, destination=m.fs.mix6.inlet_2)
m.fs.e_shell_pass7_m2 = Arc(source=m.fs.e9.outlet_2, destination=m.fs.mix6.inlet_3)
m.fs.e_shell_pass7_m3 = Arc(source=m.fs.e10.outlet_2, destination=m.fs.mix6.inlet_4)
m.fs.e_shell_pass7_m4 = Arc(source=m.fs.e11.outlet_2, destination=m.fs.mix6.inlet_5)
m.fs.e_shell_pass7_m5 = Arc(source=m.fs.e12.outlet_2, destination=m.fs.mix6.inlet_6)
m.fs.e_shell_pass7_m6 = Arc(source=m.fs.e13.outlet_2, destination=m.fs.mix6.inlet_7)

# Mixer / Separator Arc for pass 7
m.fs.e_shell_pass7_ms_Arc = Arc(source=m.fs.mix6.outlet, destination=m.fs.sep6.inlet)

# Separator Arcs for pass 7
m.fs.e_shell_pass7_s0 = Arc(source=m.fs.sep6.outlet_1, destination=m.fs.e0.inlet_2)
m.fs.e_shell_pass7_s1 = Arc(source=m.fs.sep6.outlet_2, destination=m.fs.e1.inlet_2)
m.fs.e_shell_pass7_s2 = Arc(source=m.fs.sep6.outlet_3, destination=m.fs.e2.inlet_2)
m.fs.e_shell_pass7_s3 = Arc(source=m.fs.sep6.outlet_4, destination=m.fs.e3.inlet_2)
m.fs.e_shell_pass7_s4 = Arc(source=m.fs.sep6.outlet_5, destination=m.fs.e4.inlet_2)
m.fs.e_shell_pass7_s5 = Arc(source=m.fs.sep6.outlet_6, destination=m.fs.e5.inlet_2)
m.fs.e_shell_pass7_s6 = Arc(source=m.fs.sep6.outlet_7, destination=m.fs.e6.inlet_2)

# Apply Arc constraints
pe.TransformationFactory("network.expand_arcs").apply_to(m.fs)

# Initialize Mixers and Separators
m.fs.mix0.initialize()
m.fs.mix1.initialize()
m.fs.mix2.initialize()
m.fs.mix3.initialize()
m.fs.mix4.initialize()
m.fs.mix5.initialize()
m.fs.mix6.initialize()

m.fs.sep0.initialize()
m.fs.sep1.initialize()
m.fs.sep2.initialize()
m.fs.sep3.initialize()
m.fs.sep4.initialize()
m.fs.sep5.initialize()
m.fs.sep6.initialize()


solver = pe.SolverFactory('ipopt')
solver.options = {
            "tol": 1e-5,
            "linear_solver": "ma27",
            "max_iter": 500,
        }
solver.solve(m, tee=True)


hd = 0
for idx, e in enumerate(all_elements):
    print(f'')
    print(f'-------ELEMENT #{idx}---------')
    hd += print_results_0d(e)

print('')
print(f'Total heat duty: {hd}')


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
ax.set_title('2D Model Temperature Profile')
plt.show()
