"""The LVS runset in `siphopdk/rules/siphopdk.lvs`, run through `siphopdk.lvs_check`.

A runset that only ever passes proves nothing, so each intact circuit is paired with faults that
must fail: outputs swapped, a route missing, a route moved off its pins, a metal short and an
open. The schematic is always written by PhotonForge from a component's own connections.
"""

import photonforge as pf
import pytest
from helpers import needs_klayout

import siphopdk

pytestmark = needs_klayout

RADIUS = 10.0


@pytest.fixture(scope="module")
def technology():
    tech = siphopdk.siphopdk()
    pf.config.default_technology = tech
    return tech


def lvs(component, schematic=None):
    return siphopdk.lvs_check(component, schematic=schematic, show_progress=False)


def splitter(technology, *, swap=False, drop=False, offset=0.0):
    """Grating coupler into a Y-branch: one output to a grating coupler, one to a termination."""
    c = pf.Component("splitter", technology=technology)
    gc = siphopdk.grating_coupler(technology=technology)
    term = siphopdk.termination(technology=technology)
    g_in = c.add_reference(gc).rotate(180).translate((-60, 0))
    y = c.add_reference(siphopdk.y_branch(technology=technology))
    g1 = c.add_reference(gc).translate((80, 30))
    t2 = c.add_reference(term).rotate(180).translate((80, -30))
    a, b = (t2, g1) if swap else (g1, t2)
    pairs = [(g_in["P0"], y["P0"]), (y["P1"], a["P0"]), (y["P2"], b["P0"])]
    for p, q in pairs[:2] if drop else pairs:
        c.add_reference(pf.parametric.route(port1=p, port2=q, radius=RADIUS, technology=technology))
    if offset:
        c.references[-1].translate((0, offset))
    return c


def heater(technology, *, short=False, open_=False):
    """A thermal phase shifter whose two heater terminals go to two bond pads on M1."""
    c = pf.Component("heater", technology=technology)
    ps = c.add_reference(siphopdk.thermal_phase_shifter(technology=technology))
    pad = siphopdk.bond_pad(technology=technology)
    p1 = c.add_reference(pad).translate((-40, 120))
    p2 = c.add_reference(pad).translate((140, 120))
    routes = [((ps, "HTRP"), (p1, "M1")), ((ps, "HTRN"), (p2, "M1"))]
    for t1, t2 in routes[:1] if open_ else routes:
        c.add_reference(
            pf.parametric.route_manhattan(
                terminal1=t1, terminal2=t2, width=2.0, direction1="y", technology=technology
            )
        )
    if short:
        c.add("M1", pf.Rectangle((-40, 60), (140, 62)))
    return c


def test_intact_splitter_matches(technology):
    assert lvs(splitter(technology)).verdict == "match"


@pytest.mark.parametrize(
    "fault", [dict(swap=True), dict(drop=True), dict(offset=1.0)], ids=["swapped", "missing", "off_pin"]
)
def test_optical_faults_are_caught(technology, fault):
    assert lvs(splitter(technology, **fault), schematic=splitter(technology)).verdict == "mismatch"


def test_intact_heater_matches(technology):
    """The pads' M1 and T terminals are one node, which the runset knows."""
    assert lvs(heater(technology)).verdict == "match"


@pytest.mark.parametrize("fault", [dict(short=True), dict(open_=True)], ids=["short", "open"])
def test_electrical_faults_are_caught(technology, fault):
    assert lvs(heater(technology, **fault), schematic=heater(technology)).verdict == "mismatch"


@pytest.mark.parametrize("name", siphopdk.component_names)
def test_every_cell_matches_on_its_own(name, technology):
    top = pf.Component(f"top_{name}", technology=technology)
    top.add_reference(getattr(siphopdk.component, name)(technology=technology))
    assert lvs(top).verdict == "match"


def electrode(technology, *, open_=False):
    """A travelling-wave electrode with both ends of its signal line routed to pads on M2."""
    c = pf.Component("electrode", technology=technology)
    cell = siphopdk.tw_electrode(technology=technology)
    e = c.add_reference(cell)
    pad = siphopdk.bond_pad(technology=technology)
    (x0, _), (x1, _) = cell.bounds()
    p1 = c.add_reference(pad).translate((x0 - 150, 0))
    p2 = c.add_reference(pad).translate((x1 + 150, 0))
    routes = [((e, "SIN"), (p1, "T")), ((e, "SOUT"), (p2, "T"))]
    for t1, t2 in routes[:1] if open_ else routes:
        c.add_reference(
            pf.parametric.route_manhattan(terminal1=t1, terminal2=t2, width=5.0, technology=technology)
        )
    return c


def test_electrode_with_both_ends_routed_matches(technology):
    """SIN and SOUT are one node, so both pads join one net, in the layout and the schematic."""
    assert lvs(electrode(technology)).verdict == "match"
    assert lvs(electrode(technology, open_=True), schematic=electrode(technology)).verdict == "mismatch"


def nitride(technology, *, drop=False):
    """A grating coupler routed to a silicon-to-nitride transition, butt-joined to an edge coupler."""
    c = pf.Component("nitride", technology=technology)
    gc = c.add_reference(siphopdk.grating_coupler(technology=technology)).rotate(180)
    s2n = c.add_reference(siphopdk.si_to_sin(technology=technology)).translate((60, 20))
    ec = c.add_reference(siphopdk.sin_edge_coupler(technology=technology))
    ec.connect("P0", s2n["P1"])
    if not drop:
        c.add_reference(
            pf.parametric.route(port1=gc["P0"], port2=s2n["P0"], radius=RADIUS, technology=technology)
        )
    return c


def test_nitride_path_and_butt_joint_match(technology):
    assert lvs(nitride(technology)).verdict == "match"
    assert lvs(nitride(technology, drop=True), schematic=nitride(technology)).verdict == "mismatch"


def mzi(technology, *, drop=False):
    """An MZI from two Y-branches, made with component_from_netlist, with two top-level ports."""
    waveguide = pf.AnalyticWaveguideModel(n_eff=2.59, n_group=4.29, propagation_loss=2e-4)
    y = siphopdk.component.y_branch(technology=technology)
    routes = [
        (("lower", "P1"), ("upper", "P2"), {"radius": 10, "model": waveguide}),
        (("lower", "P2"), ("upper", "P1"), {"radius": 20, "model": waveguide}),
    ]
    return pf.component_from_netlist(
        {
            "name": "MZI",
            "instances": {"lower": y, "upper": {"component": y, "origin": (0, 40)}},
            "routes": routes[:1] if drop else routes,
            "ports": [("lower", "P0"), ("upper", "P0")],
            "models": [pf.CircuitModel()],
        }
    )


def test_netlist_circuit_with_top_level_ports_matches(technology):
    """The ports a circuit exposes become pins of its top circuit, as in the schematic."""
    assert lvs(mzi(technology)).verdict == "match"
    assert lvs(mzi(technology, drop=True), schematic=mzi(technology)).verdict == "mismatch"


def test_top_level_terminals_match(technology):
    c = heater(technology)
    ps, p1, p2 = c.references[:3]
    c.add_port(ps["P0"], "IN")
    c.add_port(ps["P1"], "OUT")
    c.add_terminal(p1["T"], "HEAT_P")
    c.add_terminal(p2["T"], "HEAT_N")
    assert lvs(c).verdict == "match"
