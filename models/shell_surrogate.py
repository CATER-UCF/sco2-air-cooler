from pyomo.environ import Var, Param
from idaes.generic_models.unit_models import Feed
import math


class ShellSurrogate(Feed):
    """
    Extension of the IDAES Feed unit model. Used to calculate hconv and dP
    as functions of flow and property states.

    Work in progress...
    """
    def add_geometry(self):

        self.pi = Param(initialize=math.pi)

        # Tube parameters
        self.tube_outer_diameter = Param(initialize=1., mutable=True)
        self.st_pitch = Param(initialize=1., mutable=True)
        self.lt_pitch = Param(initialize=1., mutable=True)
        self.n_rows = Param(initialize=1., mutable=True)
        self.n_columns = Param(initialize=1., mutable=True)
        self.n_passes = Param(initialize=1., mutable=True)
        self.tube_length = Param(initialize=1., mutable=True)

        # Fin parameters
        self.fin_outer_diameter = Param(initialize=1., mutable=True)
        self.fin_pitch = Param(initialize=1., mutable=True)
        self.fin_thickness = Param(initialize=1., mutable=True)
        self.n_fins = Param(initialize=1., mutable=True)

        # Geometry expressions
        @self.Expression()
        def segment_length(b):
            return b.tube_length / b.n_passes

        @self.Expression()
        def flow_length(b):
            return b.st_pitch * b.n_passes + b.fin_outer_diameter

        @self.Expression()
        def traverse_length(b):
            return (b.lt_pitch * b.n_rows + b.fin_outer_diameter) * b.n_columns

        @self.Expression()
        def area_in(b):
            return b.segment_length * b.traverse_length

        @self.Expression()
        def volume_total(b):
            return b.area_in * b.flow_length

        @self.Expression()
        def volume_tubes(b):
            return (b.n_rows * b.n_columns) * 0.25 * b.pi * \
                   b.tube_outer_diameter ** 2 * b.tube_length

        @self.Expression()
        def volume_fins(b):
            return b.n_fins * 0.25 * b.pi * \
                   (b.fin_outer_diameter ** 2 - b.tube_outer_diameter ** 2) \
                   * b.fin_thickness

        @self.Expression()
        def volume_air(b):
            return b.volume_total - b.volume_tubes - b.volume_fins

        @self.Expression()
        def volume_porosity(b):
            return b.volume_air / b.volume_total

        @self.Expression()
        def tube_surface_area(b):
            return (b.n_rows * b.n_columns * b.tube_length - b.n_fins * b.fin_thickness) * \
                   b.pi * b.tube_outer_diameter

        @self.Expression()
        def fin_surface_area(b):
            return b.n_fins * (0.25 * b.pi * (b.fin_outer_diameter ** 2 -
                                              b.tube_outer_diameter ** 2)
                               * 2 + b.pi * b.fin_thickness *
                               b.fin_outer_diameter)

        @self.Expression()
        def wetted_surface_area(b):
            return b.tube_surface_area + b.fin_surface_area

        @self.Expression()
        def air_hydraulic_diameter(b):
            return 4 * b.volume_air / b.wetted_surface_area

    def add_common_eqs(self):
        self.v_in = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Velocity of incoming air'
        )
        self.v_max = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Maximum air velocity reached at the minimum cross sectional area'
        )
        self.mw_air = Var(
            self.flowsheet().config.time,
            initialize=1.,
            doc='Molecular weight of gas mixture in kg/mol'
        )
        self.N_Re_air = Var(
            self.flowsheet().config.time,
            initialize=100.,
            doc='Reynolds number'
        )

        @self.Constraint(
            self.flowsheet().config.time
        )
        def v_in_eq(b, t):
            return b.v_in[t] * b.area_in * \
                   b.properties[t].dens_mol_phase['Vap'] == b.properties[t].flow_mol

        @self.Constraint(
            self.flowsheet().config.time
        )
        def v_max_eq(b, t):
            return b.v_max[t] == b.v_in[t] / b.volume_porosity

        @self.Constraint(
            self.flowsheet().config.time
        )
        def mw_air_eqn(b, t):
            return b.mw_air[t] == sum(b.properties[t].mw_comp[j] *
                                      b.properties[t].mole_frac_comp[j]
                                      for j in b.properties[t].params.component_list)

        @self.Constraint(
            self.flowsheet().config.time
        )
        def N_Re_air_eqn(b, t):
            return b.N_Re_air[t] * \
                   b.properties[t].visc_d == \
                   b.air_hydraulic_diameter * b.v_max[t] * \
                   b.properties[t].dens_mol_phase['Vap'] * \
                   b.mw_air[t]
