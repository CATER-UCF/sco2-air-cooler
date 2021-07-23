from pyomo.environ import Var, Param
from idaes.generic_models.unit_models import Feed


class TubeSurrogate(Feed):
    """
    Extension of the IDAES Feed unit model. Used to calculate hconv and dP
    as functions of flow and property states.
    """
    def add_geometry(self):
        self.tube_inner_perimeter = Param(initialize=1., mutable=True)
        self.tube_length = Param(initialize=1., mutable=True)
        self.number_of_tubes = Param(initialize=1., mutable=True)
        self.tube_hydraulic_diameter = Param(initialize=1., mutable=True)
        self.tube_flow_area = Param(initialize=1., mutable=True)

    def add_hconv_eqs(self):
        self.v_tube = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Velocity'
        )
        self.N_Re_tube = Var(
            self.flowsheet().config.time,
            initialize=10000.,
            doc='Reynolds number'
        )
        self.N_Pr_tube = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Prandtl number'
        )
        self.N_Nu_tube = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Nusselt number'
        )
        self.hconv_tube = Var(
            self.flowsheet().config.time,
            initialize=1000.,
            doc='Heat transfer coefficient'
        )

        @self.Expression(
            self.flowsheet().config.time,
            doc='Saturated vapor fraction'
        )
        def tube_x(b, t):
            return b.properties[t].vapor_frac

        @self.Expression(
            self.flowsheet().config.time,
            doc='Density'
        )
        def rho(b, t):
            return b.properties[t].dens_mass

        @self.Expression(
            self.flowsheet().config.time,
            doc='Cp per mole'
        )
        def Cp_mol(b, t):
            return b.properties[t].cp_mol

        @self.Expression(
            self.flowsheet().config.time,
            doc='Dynamic viscosity'
        )
        def visc_d(b, t):
            return b.properties[t].visc_d_phase['Vap'] * b.tube_x[t] + \
                   b.properties[t].visc_d_phase['Liq'] * (1. - b.tube_x[t])

        @self.Expression(
            self.flowsheet().config.time,
            doc='Thermal conductivity'
        )
        def therm_cond(b, t):
            return b.properties[t].therm_cond_phase['Vap'] * b.tube_x[t] + \
                   b.properties[t].therm_cond_phase['Liq'] * (1. - b.tube_x[t])

        @self.Constraint(
            self.flowsheet().config.time, doc='Velocity equation'
        )
        def v_tube_eq(b, t):
            return 0.01 * b.v_tube[t] * b.number_of_tubes * \
                   b.tube_flow_area * \
                   b.rho[t] == \
                   0.01 * b.properties[t].flow_mass

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Reynolds number equation'
        )
        def N_Re_tube_eqn(b, t):
            return b.N_Re_tube[t] * \
                   b.visc_d[t] == \
                   b.tube_hydraulic_diameter * b.v_tube[t] * \
                   b.rho[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Prandtl number equation'
        )
        def N_Pr_tube_eqn(b, t):
            return b.N_Pr_tube[t] * \
                   b.therm_cond[t] * b.properties[t].mw == \
                   b.Cp_mol[t] * b.visc_d[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Nusselt number equation'
        )
        def N_Nu_tube_eq(b, t):
            return b.N_Nu_tube[t] == 0.023 * b.N_Re_tube[t] ** 0.8 * \
                   b.N_Pr_tube[t] ** 0.4

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Heat transfer coefficient equation'
        )
        def hconv_eq(b, t):
            return b.hconv_tube[t] * b.tube_hydraulic_diameter == \
                   b.N_Nu_tube[t] * \
                   b.therm_cond[t]
