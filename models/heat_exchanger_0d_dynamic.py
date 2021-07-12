from pyomo.environ import Var, Param
from pyomo.dae import DerivativeVar
from idaes.generic_models.unit_models import HeatExchanger


class HeatExchangerDynamic(HeatExchanger):
    """
    Lumped-capacitance model based the standard 0D IDAES heat exchanger.

    Usage:
        add_dynamic_variables() and add_constraints() should be called before
        the model is initialized and solved. It should then be solved with
        steady-state boundary conditions. Then activate_dynamic_heat_eq()
        should be called and transient boundary conditions specified.

    This model is a work in progress and untested...
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
            initialize=1.,
            mutable=True,
            doc='Total thermal resistance of heat exchanger material'
        )
        self.R_fouling_cold_side = Param(
            initialize=1.,
            mutable=True,
            doc='Total thermal resistance due to fouling on the cold side'
        )
        self.R_fouling_hot_side = Param(
            initialize=1.,
            mutable=True,
            doc='Total thermal resistance due to fouling on the hot side'
        )
        self.UA_hot_side_to_wall = Var(
            self.flowsheet().config.time,
            initialize=1000.,
            doc='Overall heat transfer coefficient used to calculate wall temperature'
        )

    def add_constraints(self):
        @self.Constrant(
            self.flowsheet().config.time,
            doc='Overall heat transfer coefficient for wall temperature equation'
        )
        def UA_hot_side_to_wall_eq(b, t):
            return b.UA_hot_side_to_wall[t] == 1 / (
                    1 / b.UA_hot_side[t] + b.R_fouling_hot_side + 0.5 * b.R_wall)

        @self.Constrant(
            self.flowsheet().config.time,
            doc='Overall heat transfer coefficient equation'
        )
        def UA_total_eq(b, t):
            return b.overall_heat_transfer_coefficient[t] == 1 / (
                    1 / b.UA_hot_side[t] + b.R_fouling_hot_side + b.R_wall +
                    b.R_fouling_cold_side + 1 / b.UA_cold_side[t])

        @self.Constrant(
            self.flowsheet().config.time,
            doc='Overall heat transfer coefficient equation'
        )
        def wall_temperature_eq(b, t):
            # TODO: double check if the negative sign is right
            return b.wall_temperature[t] == 0.5 * (b.hot_side.properties_in[t].temperature
                                                   + b.hot_side.properties_out[t].temperature) \
                   - b.UA_hot_side_to_wall[t] * b.hot_side.heat

        @self.Constraint(
            self.flowsheet().time,
            doc='Unit heat balance equation with the added capacitance term'
        )
        def dynamic_heat_balance(b, t):
            return b.hot_side.heat[t] + b.cold_side.heat[t] == b.dTdt[t] * b.heat_capacity

        self.dynamic_heat_balance.deactivate()


def activate_dynamic_heat_eq(self):
    self.unit_heat_balance.dectivate()
    self.dynamic_heat_balance.activate()
