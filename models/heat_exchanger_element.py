from pyomo.environ import Var, Param
from models import HeatExchangerLumpedCapacitance


class HeatExchangerElement(HeatExchangerLumpedCapacitance):

    def setup(self):
        self.tube_length = Param(initialize=1, mutable=True)
        self.internal_surface_area = Param(initialize=1, mutable=True)
        self.external_surface_area = Param(initialize=1, mutable=True)
        self.add_dynamic_variables()
        self.add_dynamic_variable_constraints()
        self.add_pressure_flow_constraints()
        self.add_surrogate_constraints()
        self.add_UA_and_flow_coefficient_constraints()

    def add_surrogate_constraints(self):

        self.hconv_tube = Var(
            self.flowsheet().config.time,
            initialize=1000.
        )

        self.dP_over_l_tube = Var(
            self.flowsheet().config.time,
            initialize=-1000.
        )

        self.hconv_shell = Var(
            self.flowsheet().config.time,
            initialize=20.
        )

        # Predictor variables
        @self.Expression(
            self.flowsheet().config.time
        )
        def tube_cp_mol(b, t):
            return 0.5 * (b.tube.properties_in[t].cp_mol +
                          b.tube.properties_out[t].cp_mol)

        @self.Expression(
            self.flowsheet().config.time
        )
        def tube_rho(b, t):
            return 0.5 * (b.tube.properties_in[t].dens_mass +
                          b.tube.properties_out[t].dens_mass)

        @self.Expression(
            self.flowsheet().config.time
        )
        def shell_temp(b, t):
            return 0.5 * (b.shell.properties_in[t].temperature +
                          b.shell.properties_out[t].temperature)

        # Surrogate models
        @self.Constraint(
            self.flowsheet().config.time
        )
        def hconv_tube_surrogate_eq(b, t):
            return b.hconv_tube[t] == -2.57411323e+03 + 1.45303372e+03 \
                   * b.tube_cp_mol[t] ** 3.19268561e-01 + \
                   -3.73594780e-02 * b.tube_cp_mol[t]

        @self.Constraint(
            self.flowsheet().config.time
        )
        def dP_tube_surrogate_eq(b, t):
            return b.dP_over_l_tube[t] == 2.19442286e+01 - \
                   2.68411214e+05 / b.tube_rho[t] + \
                   -7.17992324e-02 * b.tube_rho[t] + \
                   3.93218678e-05 * b.tube_rho[t] ** 2

        @self.Constraint(
            self.flowsheet().config.time
        )
        def hconv_shell_surrogate_eq(b, t):
            return b.hconv_shell[t] == 5.00519828 + 0.04982918 * \
                   b.shell_temp[t]

    def add_UA_and_flow_coefficient_constraints(self):

        @self.Constraint(self.flowsheet().config.time)
        def UA_hot_side_eq(b, t):
            return b.UA_hot_side[t] == b.hconv_tube[t] * b.internal_surface_area

        @self.Constraint(self.flowsheet().config.time)
        def UA_cold_side_eq(b, t):
            return b.UA_cold_side[t] == b.hconv_shell[t] * b.external_surface_area

        @self.Constraint(self.flowsheet().config.time)
        def flow_coefficient_eq(b, t):
            return b.flow_coefficient_hot_side[t] == -b.dP_over_l_tube[t] * \
                   b.tube_length / 611.6
