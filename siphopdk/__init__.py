"""SiPhoPDK: a synthetic 220 nm SOI + SiN O-band demonstration PDK for PhotonForge.

Not a foundry PDK. Every number here is either published physics, taken from an open
source, or invented for teaching. Do not tape out against it.

    import photonforge as pf
    import siphopdk

    pf.config.default_technology = siphopdk.siphopdk()
    ring = siphopdk.component.micro_ring_modulator()
"""

from . import component
from .cells import (
    bend,
    bond_pad,
    crossing,
    directional_coupler,
    grating_coupler,
    micro_ring_modulator,
    micro_ring_resonator,
    mmi_2x2,
    photodiode,
    pn_phase_shifter,
    si_to_sin,
    sin_edge_coupler,
    straight,
    strip_to_rib,
    tap,
    taper,
    termination,
    thermal_phase_shifter,
    via_stack,
    y_branch,
)
from .ladder_cells import (
    add_drop_ring,
    add_drop_ring_corner,
    contra_dc,
    grating_antenna,
    pin_attenuator,
    spiral_delay,
    tw_electrode,
)
from ._cross_section import (
    cross_section_component,
    layer_coverage,
    plot_cross_section,
    plot_stack_panels,
)
from .technology import siphopdk
from .verification import DRC_RULES, LVS_RULES, drc_check, lvs_check

__version__ = "1.0.0"

#: The technology function, under a second name.
technology = siphopdk

#: Names of all cells in ``siphopdk.component``.
component_names = tuple(sorted(n for n in dir(component) if not n.startswith("_")))

__all__ = [
    "DRC_RULES",
    "LVS_RULES",
    "__version__",
    "add_drop_ring",
    "add_drop_ring_corner",
    "bend",
    "bond_pad",
    "component",
    "component_names",
    "contra_dc",
    "cross_section_component",
    "crossing",
    "directional_coupler",
    "drc_check",
    "grating_antenna",
    "grating_coupler",
    "layer_coverage",
    "lvs_check",
    "micro_ring_modulator",
    "micro_ring_resonator",
    "mmi_2x2",
    "photodiode",
    "pin_attenuator",
    "plot_cross_section",
    "plot_stack_panels",
    "pn_phase_shifter",
    "si_to_sin",
    "sin_edge_coupler",
    "siphopdk",
    "spiral_delay",
    "straight",
    "strip_to_rib",
    "tap",
    "taper",
    "technology",
    "termination",
    "thermal_phase_shifter",
    "tw_electrode",
    "via_stack",
    "y_branch",
]
