from pyomo.environ import Var, Param
from models import HeatExchangerLumpedCapacitance


class HeatExchangerElement(HeatExchangerLumpedCapacitance):
    """
    Heat exchanger finite element for use in the two-dimensional flow network.

    Usage:
    ------
    The general procedure for this model is as follows:

        1. Add model(s) and Arcs to flowsheet
        2. Call setup()
        3. Set boundary conditions for steady-state
        4. Initialize and solve the model
        5. Add temperature disturbances and call activate_dynamic_heat_eq()
        6. Solve the model with transient conditions

    Correlations from ./surrogates are implemented with the following methods:

        hconv_tube_surrogate_eq()
        dP_tube_surrogate_eq()
        hconv_shell_surrogate_eq()
        critical_temp_in()
        critical_temp_out()

    """
    def setup(self):
        self.tube_length = Param(initialize=1, mutable=True)
        self.internal_surface_area = Param(initialize=1, mutable=True)
        self.external_surface_area = Param(initialize=1, mutable=True)
        self.critical_temp_margin = Param(initialize=-999., mutable=True)
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
            return b.hconv_tube[t] == -2.73154315e+03 + 1.23970466e+03 \
                   * b.tube_cp_mol[t] ** 3.17715435e-01 + \
                   -2.82026104e-02 * b.tube_cp_mol[t]

        @self.Constraint(
            self.flowsheet().config.time
        )
        def dP_tube_surrogate_eq(b, t):
            return b.dP_over_l_tube[t] == 8.82481388e+01 - \
                   1.81263728e+05 / b.tube_rho[t] + \
                   -2.98465078e-01 * b.tube_rho[t] + \
                   1.92815541e-04 * b.tube_rho[t] ** 2

        @self.Constraint(
            self.flowsheet().config.time
        )
        def hconv_shell_surrogate_eq(b, t):
            return b.hconv_shell[t] == 5.00519828 + 0.04982918 * \
                   b.shell_temp[t]

        @self.Expression(
            self.flowsheet().config.time
        )
        def critical_temp_in(b, t):
            return 2.55996364e+02 + \
                   7.15290022e-06 * b.tube.properties_in[t].pressure \
                   - 8.52339353e-14 * b.tube.properties_in[t].pressure ** 2

        @self.Expression(
            self.flowsheet().config.time
        )
        def critical_temp_out(b, t):
            return 2.55996364e+02 + \
                   7.15290022e-06 * b.tube.properties_out[t].pressure \
                   - 8.52339353e-14 * b.tube.properties_out[t].pressure ** 2

        # Using these constraints eliminates the AMPL evaluation errors from
        # ipopt. However, the model takes much longer to solve this way. So
        # I'm leaving this in as an option in case it's useful later. But it's
        # disabled by setting critical_temp_margin=-999.
        @self.Constraint(
            self.flowsheet().config.time,
        )
        def temp_in_above_critical(b, t):
            return b.tube.properties_in[t].temperature >= \
                   b.critical_temp_in[t] + b.critical_temp_margin

        @self.Constraint(
            self.flowsheet().config.time
        )
        def temp_out_above_critical(b, t):
            return b.tube.properties_out[t].temperature >= \
                   b.critical_temp_out[t] + b.critical_temp_margin

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
