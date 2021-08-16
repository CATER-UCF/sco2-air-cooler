"""
Flowsheet for generating a surrogate model of hconv and dP as functions of
fluid state.
"""
import numpy as np
import pyomo.environ as pe
from idaes.core import FlowsheetBlock
from idaes.power_generation.properties import FlueGasParameterBlock
from models import ShellSurrogate
import pandas as pd


def make_model(temp, press, flow_mol, dyn=True, n_pts=10):

    m = pe.ConcreteModel()
    if dyn:
        m.fs = FlowsheetBlock(default={"dynamic": True,
                                       "time_set": [0, n_pts],
                                       "time_units": pe.units.s})
    else:
        m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.prop_flue_gas = FlueGasParameterBlock()
    m.fs.feed = ShellSurrogate(default={"property_package": m.fs.prop_flue_gas})
    m.fs.feed.add_geometry()
    m.fs.feed.add_common_eqs()

    if dyn:
        m.discretizer = pe.TransformationFactory('dae.finite_difference')
        m.discretizer.apply_to(m, nfe=n_pts - 1, wrt=m.fs.time, scheme="BACKWARD")

    m.fs.feed.tube_outer_diameter = 0.0318
    m.fs.feed.st_pitch = 0.0699
    m.fs.feed.lt_pitch = 0.0605
    m.fs.feed.n_rows = 29
    m.fs.feed.n_columns = 40
    m.fs.feed.n_passes = 8
    m.fs.feed.tube_length = 195
    m.fs.feed.fin_outer_diameter = 0.0762
    m.fs.feed.fin_pitch = 0.0025
    m.fs.feed.fin_thickness = 0.0004
    m.fs.feed.n_fins = 76800 * 1160

    m.fs.feed.outlet.flow_mol_comp[:, "H2O"].fix(0.01027 * flow_mol)
    m.fs.feed.outlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * flow_mol)
    m.fs.feed.outlet.flow_mol_comp[:, "N2"].fix(0.780066026 * flow_mol)
    m.fs.feed.outlet.flow_mol_comp[:, "O2"].fix(0.209252382 * flow_mol)
    m.fs.feed.outlet.flow_mol_comp[:, "NO"].fix(0)
    m.fs.feed.outlet.flow_mol_comp[:, "SO2"].fix(0)

    m.fs.feed.outlet.pressure[:].fix(press)
    m.fs.feed.outlet.temperature[:].fix(temp)

    return m


def run_doe(npr, nt):

    # Full-factorial DOE over a range of temperatures and pressures
    t_min = 273.15
    t_max = 390
    p_min = 101325 * 0.9
    p_max = 101325 * 1.1

    temp, press = np.meshgrid(np.linspace(t_min, t_max, nt),
                              np.linspace(p_min, p_max, npr))

    temp = temp.flatten()
    press = press.flatten()

    solver = pe.SolverFactory('ipopt')
    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 500,
    }

    shell_flow = 44004.14222
    m = make_model(288.15, 101325, shell_flow, dyn=True, n_pts=nt*npr)

    for i, t in enumerate(m.fs.time):
        m.fs.feed.outlet.temperature[t].fix(temp[i])
        m.fs.feed.outlet.pressure[t].fix(press[i])

    m.fs.feed.initialize()
    solver.solve(m, tee=True)
    Re = pe.value(m.fs.feed.N_Re_air[:])
    hconv = pe.value(m.fs.feed.hconv_air[:])
    v_in = pe.value(m.fs.feed.v_in[:])
    v_max = pe.value(m.fs.feed.v_max[:])

    print(f'Re: {Re[0]}')
    print(f'hconv: {hconv[0]}')

    hd = pe.value(m.fs.feed.air_hydraulic_diameter)
    print(f'HD: {hd}')

    a_in = pe.value(m.fs.feed.area_in)
    print(f'Area in: {a_in}')

    phi = pe.value(m.fs.feed.volume_porosity)
    print(f'phi: {phi}')

    v_total = pe.value(m.fs.feed.volume_total)
    v_air = pe.value(m.fs.feed.volume_air)
    v_fins = pe.value(m.fs.feed.volume_fins)
    v_tubes = pe.value(m.fs.feed.volume_tubes)

    print(f'Total volume: {v_total}')
    print(f'Air volume: {v_air}')
    print(f'Fin volume: {v_fins}')
    print(f'Tube volume: {v_tubes}')

    l_flow = pe.value(m.fs.feed.flow_length)
    print(f'Flow length: {l_flow}')

    fin_area = pe.value(m.fs.feed.fin_surface_area)
    tube_area = pe.value(m.fs.feed.tube_surface_area)

    print(f'Fin area: {fin_area}')
    print(f'Tube area: {tube_area}')

    df = pd.DataFrame(data={
        'temperature': temp,
        'pressure': press,
        'Re': np.array(Re),
        'hconv': np.array(hconv),
        'v_in': np.array(v_in),
        'v_max': np.array(v_max)
    })
    df.to_csv(f'./data/DOE_air_p{npr}_t{nt}.csv', index=None)


if __name__ == '__main__':
    run_doe(20, 20)
