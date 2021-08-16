from pyomo.environ import Var, Param
from pyomo.dae import DerivativeVar
from idaes.generic_models.unit_models import HeatExchanger


class HeatExchangerLumpedCapacitance(HeatExchanger):
    """
    Lumped-capacitance model based the standard 0D IDAES heat exchanger.

    Usage
    -----
    The unit model should be solved in sequence like so:
        1. Add unit model to flowsheet
        2. Call add_dynamic_variables()
        3. Call add_dynamic_variable_constraints()
        4. Apply discretization to the time domain
        5. Set heat transfer parameters and apply steady-state boundary
           condition constraints
        6. Initialize and solve the model
        7. Apply transient boundary condition constraints
        8. Call activate_dynamic_heat_eq()
        9. Solve the model
    """
    def add_dynamic_variables(self):
        self.wall_temperature = Var(
            self.flowsheet().config.time,
            initialize=300.,
            doc='Average wall temperature'
        )
        self.dTdt = DerivativeVar(
            self.wall_temperature,
            wrt=self.flowsheet().config.time,
            doc='Derivative of wall temperature with respect to time'
        )
        self.heat_capacity = Param(
            initialize=1.,
            mutable=True,
            doc='Total heat capacity of wall and fin material (J/K)'
        )
        self.UA_cold_side = Var(
            self.flowsheet().config.time,
            initialize=1000.,
            doc='Overall heat transfer coefficient from the cold side'
        )
        self.UA_hot_side = Var(
            self.flowsheet().config.time,
            initialize=1000.,
            doc='Overall heat transfer coefficient from the hot side'
        )
        self.R_wall = Param(
            initialize=0.,
            mutable=True,
            doc='Total thermal resistance of heat exchanger material'
        )
        self.R_fouling_cold_side = Param(
            initialize=0.,
            mutable=True,
            doc='Total thermal resistance due to fouling on the cold side'
        )
        self.R_fouling_hot_side = Param(
            initialize=0.,
            mutable=True,
            doc='Total thermal resistance due to fouling on the hot side'
        )
        self.UA_hot_side_to_wall = Var(
            self.flowsheet().config.time,
            initialize=1000.,
            doc='Overall heat transfer coefficient used to calculate wall temperature'
        )
        self.flow_coefficient_hot_side = Var(
            self.flowsheet().config.time,
            initialize=0.,
            doc='Relates dP to mass flow'
        )
        self.flow_coefficient_cold_side = Var(
            self.flowsheet().config.time,
            initialize=0.,
            doc='Relates dP to mass flow'
        )

    def add_dynamic_variable_constraints(self):
        @self.Constraint(
            self.flowsheet().config.time,
            doc='Overall heat transfer coefficient for wall temperature equation'
        )
        def UA_hot_side_to_wall_eq(b, t):
            return b.UA_hot_side_to_wall[t] == 1 / (
                    1 / b.UA_hot_side[t] + b.R_fouling_hot_side +
                    0.5 * b.R_wall)

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Overall heat transfer coefficient equation'
        )
        def UA_total_eq(b, t):
            return b.overall_heat_transfer_coefficient[t] == 1 / (
                    1 / b.UA_hot_side[t] + b.R_fouling_hot_side + b.R_wall +
                    b.R_fouling_cold_side + 1 / b.UA_cold_side[t])

        @self.Constraint(
            self.flowsheet().config.time,
            doc='Wall temperature equation'
        )
        def wall_temperature_eq(b, t):
            return b.wall_temperature[t] == 0.5 * (
                   b.hot_side.properties_in[t].temperature +
                   b.hot_side.properties_out[t].temperature) \
                   + b.hot_side.heat[t] / b.UA_hot_side_to_wall[t]

        @self.Constraint(
            self.flowsheet().time,
            doc='Unit heat balance equation with the added capacitance term'
        )
        def dynamic_heat_balance(b, t):
            return b.hot_side.heat[t] + b.cold_side.heat[t] == \
                   -b.dTdt[t] * b.heat_capacity

        self.dynamic_heat_balance.deactivate()
        t0 = self.flowsheet().config.time.first()
        self.dTdt[:].value = 0
        self.dTdt[t0].fix(0)

    def add_pressure_flow_constraints(self):

        if self.config.hot_side_config.has_pressure_change:
            @self.Constraint(
                self.flowsheet().time,
                doc='A simple relation for flow and dP'
            )
            def flow_coefficient_hot_eq(b, t):
                return b.hot_side.deltaP[t] == \
                       -b.flow_coefficient_hot_side[t] * \
                       b.hot_side.properties_in[t].flow_mass

        if self.config.cold_side_config.has_pressure_change:
            @self.Constraint(
                self.flowsheet().time,
                doc='A simple relation for flow and dP'
            )
            def flow_coefficient_cold_eq(b, t):
                return b.cold_side.deltaP[t] == \
                       -b.flow_coefficient_cold_side[t] * \
                       b.cold_side.properties_in[t].flow_mass

    def activate_dynamic_heat_eq(self):
        self.unit_heat_balance.deactivate()
        self.dynamic_heat_balance.activate()
