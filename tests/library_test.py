"""Every cell builds, carries the interface it promises, and answers a solver offline.

These are the checks a PDK owes its users, written so they fail rather than print. The old
`test_cells.py` printed a table for a human to read, which catches nothing in CI.
"""

import numpy as np
import photonforge as pf
import pytest

import siphopdk
from siphopdk import cells
from siphopdk._recognition import has_pin_marker

# Every cell and the interface it is contracted to expose: ports that take a pin, then
# terminals, then ports that do not take a pin. A cell that gains or loses one of these has
# changed its interface, which is the thing a PDK must never do quietly - and the split between
# the first and third groups is the rule the LVS pin markers follow, so a port that changes
# classification shows up here too.
CATALOG = {
    "straight": (("P0", "P1"), (), ()),
    "bend": (("P0", "P1"), (), ()),
    "taper": (("P0", "P1"), (), ()),
    "strip_to_rib": (("P0", "P1"), (), ()),
    "crossing": (("P0", "P1", "P2", "P3"), (), ()),
    "y_branch": (("P0", "P1", "P2"), (), ()),
    "mmi_2x2": (("P0", "P1", "P2", "P3"), (), ()),
    "directional_coupler": (("P0", "P1", "P2", "P3"), (), ()),
    "grating_coupler": (("P0",), (), ("FIBER",)),
    "tap": (("P0", "P1", "P2", "P3"), (), ()),
    "termination": (("P0",), (), ()),
    "photodiode": (("P0",), ("ANODE", "CATHODE"), ("OUT",)),
    "thermal_phase_shifter": (("P0", "P1"), ("HTRN", "HTRP"), ()),
    "pn_phase_shifter": (("P0", "P1"), ("ANODE", "CATHODE"), ()),
    "micro_ring_modulator": (("P0", "P1"), ("ANODE", "CATHODE", "HTRN", "HTRP"), ("DRIVE",)),
    "via_stack": ((), ("M1", "M2"), ()),
    "bond_pad": ((), ("M1", "T"), ()),
    "si_to_sin": (("P0", "P1"), (), ()),
    "sin_edge_coupler": (("P0",), (), ("FACET",)),
    "micro_ring_resonator": (("P0", "P1"), ("HTRN", "HTRP"), ()),
}


@pytest.fixture(scope="module", autouse=True)
def technology():
    tech = siphopdk.technology()
    pf.config.default_technology = tech
    return tech


@pytest.fixture(scope="module")
def frequencies():
    return pf.C_0 / np.linspace(1.29, 1.33, 41)


@pytest.mark.parametrize("name", sorted(CATALOG))
def test_cell_builds_with_its_declared_interface(name):
    component = getattr(cells, name)()
    optical, terminals, electrical = CATALOG[name]
    ports = set(component.ports)
    assert set(optical) | set(electrical) == ports, (
        f"{name} exposes {sorted(ports)}, the catalog says "
        f"{sorted(set(optical) | set(electrical))}"
    )
    assert set(terminals) == set(component.terminals)


@pytest.mark.parametrize("name", sorted(CATALOG))
def test_every_port_is_classified_the_way_the_pin_rule_expects(name):
    """An electrical port gets no pin marker; an optical one does. Both halves matter.

    The LVS runset reads these markers: a pin on an electrical port, or a missing pin on an
    optical one, would change the extracted netlist.
    """
    component = getattr(cells, name)()
    optical, _, electrical = CATALOG[name]
    for port_name in optical:
        assert has_pin_marker(component.ports[port_name]), f"{name}.{port_name} should get a pin"
    for port_name in electrical:
        assert not has_pin_marker(component.ports[port_name]), (
            f"{name}.{port_name} is electrical or a fibre facet and must not get a pin"
        )


@pytest.mark.parametrize("name", sorted(n for n in CATALOG if CATALOG[n][0]))
def test_every_optical_cell_answers_a_solver_offline(name, frequencies):
    """No cloud, no FDTD: a PDK cell that cannot return an S matrix on a laptop is not usable.

    Passivity is checked as a sanity bound rather than an identity, because a tap, a
    termination and a photodiode are all meant to lose light.
    """
    component = getattr(cells, name)()
    s = component.s_matrix(frequencies)
    assert s.elements, f"{name} returned an empty S matrix"
    for key, value in s.elements.items():
        assert np.all(np.isfinite(value)), f"{name} {key} has non-finite entries"
        assert np.all(np.abs(value) <= 1.0 + 1e-6), f"{name} {key} exceeds unity"


def test_ring_channels_land_on_the_declared_grid():
    """The radii and the channel wavelengths are one plan, and must agree.

    `CHANNEL_RADII` is what gets drawn and `CHANNEL_WAVELENGTHS` is what the simulation and the
    documentation quote. They are derived from the same spacing, so this checks the derivation
    rather than a table: a ring built at radius k has to resonate at wavelength k.
    """
    from siphopdk import _calibration as cal

    lam = np.linspace(1.3085, 1.3165, 1601)
    freqs = pf.C_0 / lam
    for index, (radius, target) in enumerate(zip(cal.CHANNEL_RADII, cal.CHANNEL_WAVELENGTHS)):
        ring = cells.micro_ring_modulator(ring_radius=radius)
        through = np.abs(ring.s_matrix(freqs)[("P0@0", "P1@0")]) ** 2
        found = lam[np.argmin(np.where(np.abs(lam - target) < 0.0004, through, np.inf))]
        assert abs(found - target) * 1e3 < 0.05, (
            f"channel {index + 1} resonates at {found * 1e3:.3f} nm, "
            f"the plan says {target * 1e3:.3f} nm"
        )


@pytest.mark.parametrize("port_spec,n_eff", [("SI_STRIP_O", 2.5932), ("SI_RIB_O", 2.7856),
                                             ("SIN_STRIP_O", 1.6429)])
def test_a_guide_is_modelled_with_its_own_index(port_spec, n_eff):
    """A cell drawn on one layer must not be modelled as another.

    `straight(port_spec="SIN_STRIP_O")` draws on nitride and used to be handed silicon's
    index, because the model helper took a two-state `rib=` flag and the PDK has three guides.
    Nothing anywhere reported it: the layout was right, the DRC was clean, and the delay was
    out by seventy percent.
    """
    import re

    component = cells.straight(length=100.0, port_spec=port_spec)
    found = float(re.search(r"n_eff=([-0-9.]+)", repr(component.active_model)).group(1))
    assert abs(found - n_eff) < 1e-3, f"{port_spec} modelled with n_eff {found}"


def test_technology_layers_and_ports_are_stable():
    """The layer table and the port specs are the PDK's ABI. A change here breaks every GDS."""
    tech = siphopdk.technology()
    assert len(tech.layers) == 30
    for name, number in (("SI", 1), ("SLAB", 2), ("GRAT", 3), ("GE", 8), ("HTR", 30),
                         ("M1", 50), ("M2", 52), ("PAD", 60),
                         ("DEVREC", 100), ("PINREC", 101), ("PINREC_E", 102)):
        assert tuple(tech.layers[name].layer) == (number, 0), f"{name} moved"
    for spec in ("SI_STRIP_O", "SI_RIB_O"):
        assert spec in tech.ports
