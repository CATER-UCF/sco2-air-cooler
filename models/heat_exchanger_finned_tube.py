from models import HeatExchangerLumpedCapacitance
from pyomo.environ import Var, Param


class HeatExchangerFinnedTube(HeatExchangerLumpedCapacitance):

    def add_geometry(self):
        self.tube_inner_perimeter = Param(initialize=1., mutable=True)
        self.tube_length = Param(initialize=1., mutable=True)
        self.number_of_tubes = Param(initialize=1., mutable=True)
        self.tube_hydraulic_diameter = Param(initialize=1., mutable=True)
        self.tube_flow_area = Param(initialize=1., mutable=True)

    def add_hconv_eqs(self):

        self.v_tube_in = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Velocity at tube inlet'
        )
        self.v_tube_out = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Velocity at tube exit'
        )
        self.N_Re_tube_in = Var(
            self.flowsheet().config.time,
            initialize=10000.,
            doc='Reynolds number at tube side inlet'
        )
        self.N_Re_tube_out = Var(
            self.flowsheet().config.time,
            initialize=10000.,
            doc='Reynolds number at tube side exit'
        )
        self.N_Pr_tube_in = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Prandtl number at the tube side inlet'
        )
        self.N_Pr_tube_out = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Prandtl number at the tube side exit'
        )
        self.N_Nu_tube_in = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Nusselt number at the tube side inlet'
        )
        self.N_Nu_tube_out = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Nusselt number at the tube side exit'
        )
        self.H_tube_in = Var(
            self.flowsheet().config.time,
            initialize=1000.,
            doc='HTC at tube inlet'
        )
        self.H_tube_out = Var(
            self.flowsheet().config.time,
            initialize=1000.,
            doc='HTC at tube exit'
        )

        @self.Expression(
            self.flowsheet().config.time,
            doc='Saturated vapor fraction at tube inlet'
        )
        def x_in(b, t):
            return b.tube.properties_in[t].vapor_frac

        @self.Expression(
            self.flowsheet().config.time,
            doc='Saturated vapor fraction at tube exit'
        )
        def x_out(b, t):
            return b.tube.properties_out[t].vapor_frac

        @self.Expression(
            self.flowsheet().config.time,
            doc='Density at tube inlet'
        )
        def rho_in(b, t):
            return b.tube.properties_in[t].dens_mass

        @self.Expression(
            self.flowsheet().config.time,
            doc='Density at tube exit'
        )
        def rho_out(b, t):
            return b.tube.properties_out[t].dens_mass

        @self.Expression(
            self.flowsheet().config.time,
            doc='Cp per mole at tube inlet'
        )
        def Cp_mol_in(b, t):
            return b.tube.properties_in[t].cp_mol

        @self.Expression(
            self.flowsheet().config.time,
            doc='Cp per mole at tube exit'
        )
        def Cp_mol_out(b, t):
            return b.tube.properties_out[t].cp_mol

        @self.Expression(
            self.flowsheet().config.time,
            doc='Dynamic viscosity at tube inlet'
        )
        def visc_d_in(b, t):
            return b.tube.properties_in[t].visc_d_phase['Vap'] * b.x_in[t] + \
                   b.tube.properties_in[t].visc_d_phase['Liq'] * (1. - b.x_in[t])

        @self.Expression(
            self.flowsheet().config.time,
            doc='Dynamic viscosity at tube exit'
        )
        def visc_d_out(b, t):
            return b.tube.properties_out[t].visc_d_phase['Vap'] * b.x_out[t] + \
                   b.tube.properties_out[t].visc_d_phase['Liq'] * (1. - b.x_out[t])

        @self.Expression(
            self.flowsheet().config.time,
            doc='Thermal conductivity at tube inlet'
        )
        def therm_cond_in(b, t):
            return b.tube.properties_in[t].therm_cond_phase['Vap'] * b.x_in[t] + \
                   b.tube.properties_in[t].therm_cond_phase['Liq'] * (1. - b.x_in[t])

        @self.Expression(
            self.flowsheet().config.time,
            doc='Thermal conductivity at tube exit'
        )
        def therm_cond_out(b, t):
            return b.tube.properties_out[t].therm_cond_phase['Vap'] * b.x_out[t] + \
                   b.tube.properties_out[t].therm_cond_phase['Liq'] * (1. - b.x_out[t])

        @self.Constraint(
            self.flowsheet().config.time, doc='Velocity at tube side inlet'
        )
        def v_tube_in_eq(b, t):
            return 0.01 * b.v_tube_in[t] * b.number_of_tubes * \
                   b.tube_flow_area * \
                   b.rho_in[t] == \
                   0.01 * b.tube.properties_in[t].flow_mass

        @self.Constraint(
            self.flowsheet().config.time, doc='Velocity at tube side exit'
        )
        def v_tube_out_eq(b, t):
            return 0.01 * b.v_tube_out[t] * b.number_of_tubes * \
                   b.tube_flow_area * \
                   b.rho_out[t] == \
                   0.01 * b.tube.properties_out[t].flow_mass

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Reynolds number equation for tube side inlet'
        )
        def N_Re_tube_in_eqn(b, t):
            return b.N_Re_tube_in[t] * \
                   b.visc_d_in[t] == \
                   b.tube_hydraulic_diameter * b.v_tube_in[t] * \
                   b.rho_in[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Reynolds number equation for tube side exit'
        )
        def N_Re_tube_out_eqn(b, t):
            return b.N_Re_tube_out[t] * \
                   b.visc_d_out[t] == \
                   b.tube_hydraulic_diameter * b.v_tube_out[t] * \
                   b.rho_out[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Prandtl number equation for tube side inlet'
        )
        def N_Pr_tube_in_eqn(b, t):
            return b.N_Pr_tube_in[t] * \
                   b.therm_cond_in[t] * b.tube.properties_in[t].mw == \
                   b.Cp_mol_in[t] * b.visc_d_in[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Prandtl number equation for tube side exit'
        )
        def N_Pr_tube_out_eqn(b, t):
            return b.N_Pr_tube_out[t] * \
                   b.therm_cond_out[t] * b.tube.properties_out[t].mw == \
                   b.Cp_mol_out[t] * b.visc_d_out[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Nusselt number at tube side inlet'
        )
        def N_Nu_tube_in_eq(b, t):
            return b.N_Nu_tube_in[t] == 0.023 * b.N_Re_tube_in[t] ** 0.8 * \
                   b.N_Pr_tube_in[t] ** 0.4

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Nusselt number at tube side exit'
        )
        def N_Nu_tube_out_eq(b, t):
            return b.N_Nu_tube_out[t] == 0.023 * b.N_Re_tube_out[t] ** 0.8 * \
                   b.N_Pr_tube_out[t] ** 0.4

        @self.Constraint(
            self.flowsheet().config.time,
            doc='HTC at tube side inlet'
        )
        def H_tube_in_eq(b, t):
            return b.H_tube_in[t] * b.tube_hydraulic_diameter == \
                   b.N_Nu_tube_in[t] * \
                   b.therm_cond_in[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='HTC at tube side exit'
        )
        def H_tube_out_eq(b, t):
            return b.H_tube_out[t] * b.tube_hydraulic_diameter == \
                   b.N_Nu_tube_out[t] * \
                   b.therm_cond_out[t]

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Tube UA calculated from the average inlet and exit HTC'
        )
        def UA_hot_side_eq(b, t):
            return b.UA_hot_side[t] == (b.H_tube_in[t] +
                                    b.H_tube_out[t]) * \
                   b.tube_inner_perimeter * \
                   b.tube_length * b.number_of_tubes

    def add_dP_eqs(self):

        self.f_tube_in = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Darcy friction factor at tube inlet'
        )
        self.f_tube_out = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Darcy friction factor at tube exit'
        )
        self.dP_tube_in = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Pressure loss calculated using fluid properties at the tube inlet'
        )
        self.dP_tube_out = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Pressure loss calculated using fluid properties at the tube exit'
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Friction factor at tube inlet using the Blasius correlation'
        )
        def friction_factor_tube_in_eq(b, t):
            return b.f_tube_in[t] * b.N_Re_tube_in[t] ** 0.25 == 0.3164

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Friction factor at tube exit using the Blasius correlation'
        )
        def friction_factor_tube_out_eq(b, t):
            return b.f_tube_out[t] * b.N_Re_tube_out[t] ** 0.25 == 0.3164

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Darcy-Weisbach equation using fluid properties at the tube inlet'
        )
        def dP_tube_in_eq(b, t):
            return b.dP_tube_in[t] == -0.5 * b.f_tube_in[t] * b.rho_in[t] * \
                                       b.v_tube_in[t] ** 2 * b.tube_length /\
                                       b.tube_hydraulic_diameter

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Darcy-Weisbach equation using fluid properties at the tube exit'
        )
        def dP_tube_out_eq(b, t):
            return b.dP_tube_out[t] == -0.5 * b.f_tube_out[t] * b.rho_out[t] * \
                                       b.v_tube_out[t] ** 2 * b.tube_length /\
                                       b.tube_hydraulic_diameter

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Tube pressure loss constraint'
        )
        def tube_pressure_loss_eq(b, t):
            return b.tube.properties_out[t].pressure == \
                   b.tube.properties_in[t].pressure \
                    + b.dP_tube_in[t]
