from pyomo.environ import Var, Param, log
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
        self.pipe_roughness = Param(initialize=0.00002, mutable=True)

    def add_common_eqs(self):
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
        self.f_tube = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Darcy friction factor'
        )
        self.dP_over_l = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Pressure loss per length'
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
            doc='Darcy-Weisbach equation'
        )
        def dP_tube_eq(b, t):
            return b.dP_over_l[t] == -0.5 * b.f_tube[t] * b.rho[t] * \
                                     b.v_tube[t] ** 2 /\
                                     b.tube_hydraulic_diameter

    def use_churchill(self):
        @self.Constraint(
            self.flowsheet().config.time,
            doc='Friction factor using the Churchill correlation'
        )
        def friction_factor_eq(b, t):
            return b.f_tube[t] == \
                   8 * ((8.0 / b.N_Re_tube[t]) ** 12 + \
                        1.0 / ((2.456 * log(
                        1.0 / ((7.0 / b.N_Re_tube[t]) ** 0.9 + 0.27 * b.pipe_roughness / b.tube_hydraulic_diameter))) ** 16 \
                               + (37530 / b.N_Re_tube[t]) ** 16) ** 1.5) ** (1 / 12)

    def use_blasius(self):

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Friction factor using the Blasius correlation'
        )
        def friction_factor_eq(b, t):
            return b.f_tube[t] * b.N_Re_tube[t] ** 0.25 == 0.3164

    def use_dittus_boelter(self):
        self.add_common_eqs()

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Dittus-Boelter correlation'
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

    def use_gnielinski(self):

        self.add_common_eqs()

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Gneilinski correlation'
        )
        def N_Nu_tube_eq(b, t):
            return b.N_Nu_tube[t] == (b.f_tube[t] / 8) * \
                                     (b.N_Re_tube[t] - 1000) * b.N_Pr_tube[t] \
                                     / (1 + 12.7 * (b.f_tube[t] / 8) ** 0.5 * (
                                     b.N_Pr_tube[t] ** (2/3) - 1))

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Heat transfer coefficient equation'
        )
        def hconv_eq(b, t):
            return b.hconv_tube[t] * b.tube_hydraulic_diameter == \
                   b.N_Nu_tube[t] * \
                   b.therm_cond[t]
