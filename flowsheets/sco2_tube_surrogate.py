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
import pandas as pd


def make_model(sco2_t, sco2_p, sco2_mol_flow, dyn=True, n_pts=10):

    m = pe.ConcreteModel()
    if dyn:
        m.fs = FlowsheetBlock(default={"dynamic": True,
                                       "time_set": [0, n_pts],
                                       "time_units": pe.units.s})
    else:
        m.fs = FlowsheetBlock(default={"dynamic": False})

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


def run_doe(npr, nt, include_low_temperatures=False):

    # Full-factorial DOE over a range of temperatures and pressures
    t_min = 300
    t_max = 390
    p_min = 7450000
    p_max = 7760000
    p_space = np.linspace(p_min, p_max, npr)
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
        print(f'Solution found... P={p}, T*={t_star[i]}')

    """
    Now setup our DOE and make sure each T*, P pair is included. 
    """

    n_cheb1 = nt // 4
    cheb1 = np.polynomial.chebyshev.chebpts1(n_cheb1)
    cheb1 = cheb1 / cheb1[-1]

    cheb2 = np.polynomial.chebyshev.chebpts1(nt - n_cheb1 + 1)
    cheb2 = cheb2 / cheb2[-1]

    presss = np.zeros(npr * nt)
    temps = np.zeros(npr * nt)
    hconvs = np.zeros(npr * nt)
    cp_mols = np.zeros(npr * nt)
    dps = np.zeros(npr * nt)
    rhos = np.zeros(npr * nt)
    Res = np.zeros(npr * nt)
    fricts = np.zeros(npr * nt)

    solver = pe.SolverFactory('ipopt')
    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 500,
    }

    for i, p in enumerate(p_space):

        h_min = swco2.htpx(T=t_min * pe.units.K, P=p * pe.units.Pa)
        h_max = swco2.htpx(T=t_max * pe.units.K, P=p * pe.units.Pa)
        h_star = swco2.htpx(T=t_star[i] * pe.units.K, P=p * pe.units.Pa)

        h_range1 = h_min + (cheb1 + 1) / 2 * (h_star - h_min)
        h_range2 = h_star + (cheb2 + 1) / 2 * (h_max - h_star)
        h_full = np.concatenate((h_range1, h_range2[1:]))

        i_start = i * nt
        i_end = (i + 1) * nt

        m = make_model(384.35, p, 13896.84163, dyn=True, n_pts=nt)

        for j, t in enumerate(m.fs.time):
            m.fs.feed.outlet.enth_mol[t].fix(h_full[j])

        m.fs.feed.initialize()
        solver.solve(m, tee=True)

        temp = np.array(pe.value(m.fs.feed.properties[:].temperature))
        hconv = np.array(pe.value(m.fs.feed.hconv_tube[:]))
        dp = np.array(pe.value(m.fs.feed.dP_over_l[:]))
        cp_mol = np.array(pe.value(m.fs.feed.properties[:].cp_mol))
        rho = np.array(pe.value(m.fs.feed.properties[:].dens_mass))
        Re = np.array(pe.value(m.fs.feed.N_Re_tube[:]))
        frict = np.array(pe.value(m.fs.feed.f_tube[:]))

        presss[i_start: i_end] = np.ones_like(h_full) * p
        temps[i_start: i_end] = temp
        hconvs[i_start: i_end] = hconv
        dps[i_start: i_end] = dp
        cp_mols[i_start: i_end] = cp_mol
        rhos[i_start: i_end] = rho
        Res[i_start: i_end] = Re
        fricts[i_start: i_end] = frict

    fig, ax = plt.subplots(2)
    ax[0].plot(temps, hconvs, '.')
    ax[1].plot(presss, hconvs, '.')
    ax[0].set_xlabel('Temperature (K)')
    ax[1].set_xlabel('Pressure (Pa)')

    fig, ax = plt.subplots(2)
    ax[0].plot(temps, dps, '.')
    ax[1].plot(presss, dps, '.')
    ax[0].set_xlabel('Temperature (K)')
    ax[1].set_xlabel('Pressure (Pa)')

    df = pd.DataFrame(data={
        'temperature': temps,
        'pressure': presss,
        'hconv': hconvs,
        'dP_over_l': dps,
        'cp_mol': cp_mols,
        'rho': rhos,
        'Re': Res,
        'friction_factor': fricts
    })
    df.to_csv(f'./data/DOE_p{npr}_t{nt}.csv', index=None)

    df_star = pd.DataFrame(data={
        'pressure': p_space,
        'T*': t_star
    })
    df_star.to_csv(f'./data/p{npr}_t_star.csv', index=None)

    plt.show()


if __name__ == '__main__':
    run_doe(30, 99)
