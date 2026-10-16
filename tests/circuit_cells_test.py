"""The cells in `siphopdk/ladder_cells.py`: ports, names, offline models, DRC, exact spiral length."""

import warnings

import numpy as np
import photonforge as pf
import pytest

import siphopdk
from siphopdk import ladder_cells as lc

from helpers import drc_counts, needs_klayout


@pytest.fixture(scope="module", autouse=True)
def technology():
    tech = siphopdk.technology()
    pf.config.default_technology = tech
    return tech


RING_CASES = {
    "default": dict(),
    "small_passive": dict(ring_radius=3.0, gap_through=0.15, gap_drop=0.15, heater=False),
    "small_heated": dict(ring_radius=3.0, gap_through=0.15, gap_drop=0.15),
    "asymmetric": dict(ring_radius=12.0, gap_through=0.18, gap_drop=0.22),
    "large": dict(ring_radius=15.0, gap_through=0.16, gap_drop=0.16),
    "leads_down": dict(ring_radius=5.0, lead_direction=270.0),
    "leads_sideways": dict(ring_radius=8.0, lead_direction=180.0),
}


def on_the_edge_and_pointing_in(cell):
    """Every optical port sits on the cell's edge with its input direction pointing into the cell."""
    (x0, y0), (x1, y1) = cell.bounds()
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    for name, port in cell.ports.items():
        x, y = port.center
        d = np.array([np.cos(np.radians(port.input_direction)), np.sin(np.radians(port.input_direction))])
        inward = np.dot(np.array([cx - x, cy - y]), d) > 0
        on_edge = min(abs(x - x0), abs(x - x1), abs(y - y0), abs(y - y1)) < 0.15   # 0.2 um pin markers
        assert inward and on_edge, f"{cell.name}.{name} at ({x:.2f}, {y:.2f}) is not an outward-facing edge port"


def drc_clean(cell, extra=None):
    """{rule: count} of the SiPhoPDK deck on the cell, with any extra references beside it."""
    return drc_counts(cell, *(extra or []))


# ---------------------------------------------------------------------------------------
# add_drop_ring
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("case", sorted(RING_CASES))
def test_add_drop_ring_interface(case):
    cell = lc.add_drop_ring(**RING_CASES[case])
    assert sorted(cell.ports) == ["P0", "P1", "P2", "P3"]
    assert sorted(cell.terminals) == ([] if RING_CASES[case].get("heater") is False else ["HTRN", "HTRP"])
    on_the_edge_and_pointing_in(cell)


@needs_klayout
@pytest.mark.parametrize("case", sorted(RING_CASES))
def test_add_drop_ring_is_drc_clean(case):
    assert drc_clean(lc.add_drop_ring(**RING_CASES[case])) == {}


def test_add_drop_ring_names_follow_every_parameter():
    """Two cells that differ in any parameter must differ in name, and the old names must not move."""
    variants = [dict(), dict(bus_extension=20.0), dict(lead_direction=180.0), dict(lead_spread=120.0),
                dict(heater_width=1.0), dict(temperature=350.0), dict(port_spec="SI_WIDE_O"),
                dict(heater=False), dict(ring_radius=6.0), dict(gap_through=0.2), dict(gap_drop=0.2)]
    names = [lc.add_drop_ring(**kw).name for kw in variants]
    assert len(set(names)) == len(names), names
    assert names[0] == "adr_r5p0_gt0p18_gd0p18_h1"


@needs_klayout
@pytest.mark.parametrize("case", ["default", "small_passive", "asymmetric"])
def test_add_drop_ring_devrec_stops_at_the_ports(case):
    """DEVREC ends at the port plane, so butt-joined rings and straights do not overlap (DEV.1)."""
    kw = RING_CASES[case]
    cell = lc.add_drop_ring(**kw)
    (dx0, _), (dx1, _) = cell.get_structures()[(100, 0)][0].bounds()
    half = max(abs(p.center[0]) for p in cell.ports.values())
    assert dx0 == pytest.approx(-half, abs=1e-6) and dx1 == pytest.approx(half, abs=1e-6)
    top = pf.Component("t", siphopdk.technology())
    a = top.add_reference(cell)
    b = top.add_reference(lc.add_drop_ring(**kw))
    b.connect("P0", a["P1"])
    s = top.add_reference(siphopdk.straight(length=10.0))
    s.connect("P0", b["P1"])
    assert drc_counts(top) == {}


def test_add_drop_ring_refuses_a_gap_the_deck_forbids():
    with pytest.raises(ValueError, match="SI.S.1"):
        lc.add_drop_ring(ring_radius=3.0, gap_through=0.07, gap_drop=0.07)


@pytest.mark.parametrize("radius", [3.0, 5.0, 10.0])
def test_add_drop_ring_drop_peaks_where_the_through_port_dips(radius):
    cell = lc.add_drop_ring(ring_radius=radius, gap_through=0.16, gap_drop=0.16, heater=False)
    lam = np.linspace(1.30, 1.32, 4001)
    s = cell.s_matrix(pf.C_0 / lam)
    through = np.abs(s[("P0@0", "P1@0")]) ** 2
    drop = np.abs(s[("P0@0", "P2@0")]) ** 2
    assert abs(lam[np.argmin(through)] - lam[np.argmax(drop)]) < 2e-5     # 20 pm
    assert 10 * np.log10(through.min()) < -15 and drop.max() > 0.5
    assert np.abs(s[("P0@0", "P3@0")]).max() < 1e-9                      # nothing leaks to the add port
    assert (np.abs(s[("P0@0", "P1@0")]) ** 2 + drop).max() <= 1 + 1e-9


def test_add_drop_ring_temperature_moves_the_resonance():
    lam = np.linspace(1.30, 1.32, 4001)
    peak = []
    for t in (293.0, 303.0):
        cell = lc.add_drop_ring(ring_radius=5.0, gap_through=0.16, gap_drop=0.16, heater=False, temperature=t)
        peak.append(lam[np.argmax(np.abs(cell.s_matrix(pf.C_0 / lam)[("P0@0", "P2@0")]) ** 2)])
    assert 0.4e-3 < peak[1] - peak[0] < 0.8e-3        # about 575 pm for 10 K at dn/dT = 1.86e-4


@pytest.mark.parametrize("radius", [4.0, 5.0])
def test_add_drop_ring_pads_clear_the_far_bus(radius):
    """The heater landing pad sits beyond the bus the lead crossed, not on top of it."""
    cell = lc.add_drop_ring(ring_radius=radius, gap_through=0.2, gap_drop=0.2)
    y_drop = radius + 0.4 + 0.2
    for name in ("HTRP", "HTRN"):
        cx, cy = cell.terminals[name].center()
        assert cy - 2.6 >= y_drop + 0.2 + 0.5 - 1e-6, f"{name} pad reaches down to {cy - 2.6:.2f} um, the bus is at {y_drop:.2f}"


# ---------------------------------------------------------------------------------------
# add_drop_ring_corner, beside a crossing
# ---------------------------------------------------------------------------------------
def test_corner_cell_interface_and_netlist_with_a_crossing():
    cell = lc.add_drop_ring_corner(ring_radius=10.0)
    assert sorted(cell.ports) == ["P0", "P1", "P2", "P3"] and sorted(cell.terminals) == ["HTRN", "HTRP"]
    top = pf.Component("t", siphopdk.technology())
    top.add_reference(siphopdk.crossing())          # instance 0, at the origin
    top.add_reference(pf.Reference(cell))           # instance 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")             # the unconnected outer ports are reported as a warning
        connections = top.get_netlist()["connections"]
    found = {frozenset({(a[0], a[1]), (b[0], b[1])}) for a, b in connections}
    # the crossing's right arm meets the row start, its top arm the column start
    assert found == {frozenset({(0, "P2"), (1, "P0")}), frozenset({(0, "P3"), (1, "P2")})}


@needs_klayout
@pytest.mark.parametrize("radius", [6.0, 8.0, 10.0, 12.0])
def test_corner_cell_beside_a_crossing_is_drc_clean(radius):
    cell = lc.add_drop_ring_corner(ring_radius=radius)
    assert drc_clean(cell, [pf.Reference(siphopdk.crossing())]) == {}


def test_corner_cell_refuses_a_ring_too_close_to_the_crossing():
    with pytest.raises(ValueError):
        lc.add_drop_ring_corner(ring_radius=5.0)


def test_corner_cell_drops_from_the_row_into_the_column():
    cell = lc.add_drop_ring_corner(ring_radius=10.0, gap_row=0.16, gap_column=0.16, heater=False)
    lam = np.linspace(1.30, 1.32, 4001)
    s = cell.s_matrix(pf.C_0 / lam)
    through = np.abs(s[("P0@0", "P1@0")]) ** 2           # along the row
    drop = np.abs(s[("P0@0", "P2@0")]) ** 2              # row start -> column start
    assert abs(lam[np.argmin(through)] - lam[np.argmax(drop)]) < 2e-5
    assert drop.max() > 0.5 and np.abs(s[("P0@0", "P3@0")]).max() < 1e-9


# ---------------------------------------------------------------------------------------
# spiral_delay
# ---------------------------------------------------------------------------------------
SPIRAL_CASES = [dict(length=4400.0), dict(length=400.0, separation=5.0), dict(length=12345.678),
                dict(length=20000.0, separation=6.0, port_spec="SIN_STRIP_O", radius=60.0)]


@pytest.mark.parametrize("kw", SPIRAL_CASES, ids=lambda kw: f"L{kw['length']:g}")
def test_spiral_has_the_length_it_was_asked_for(kw):
    cell = lc.spiral_delay(**kw)
    assert abs(cell.properties["path_length_um"] - kw["length"]) < 2e-3        # 2 nm
    assert abs(float(pf.route_length(cell)) - kw["length"]) < 2e-3
    assert sorted(cell.ports) == ["P0", "P1"]
    p0, p1 = cell.ports["P0"], cell.ports["P1"]
    assert p0.input_direction == p1.input_direction == 0.0 and p0.center[0] == p1.center[0] == 0.0
    assert abs(p0.center[1] - p1.center[1]) == pytest.approx(kw.get("separation", 4.0))


@needs_klayout
@pytest.mark.parametrize("kw", SPIRAL_CASES, ids=lambda kw: f"L{kw['length']:g}")
def test_spiral_is_drc_clean(kw):
    assert drc_clean(lc.spiral_delay(**kw)) == {}


def test_spiral_model_loses_what_its_length_says():
    cell = lc.spiral_delay(length=4400.0)
    t = np.abs(cell.s_matrix(pf.C_0 / np.array([1.31]))[("P0@0", "P1@0")]) ** 2
    assert 0.6 < float(t[0]) < 0.9                        # 4.4 mm at 2 dB/cm plus the bends, about -1.1 dB


def test_spiral_names_follow_every_parameter():
    variants = [dict(length=1000.0), dict(length=1000.0, euler_fraction=0.3), dict(length=1000.0, lead=2.0),
                dict(length=1000.0, separation=5.0), dict(length=1000.0, radius=8.0)]
    names = [lc.spiral_delay(**kw).name for kw in variants]
    assert len(set(names)) == len(names), names


# =======================================================================================
# pin_attenuator, grating_antenna, tw_electrode, contra_dc
# =======================================================================================
from siphopdk import _calibration as cal


def test_strip_index_table_is_on_the_calibration_footing():
    """The width table's 400 nm row is the calibration table's strip, and the dispersion reproduces the solved 1290 / 1330 nm."""
    n, g = lc.strip_index(0.40)
    assert float(n) == pytest.approx(cal.N_EFF_STRIP, abs=1e-4) and g == pytest.approx(cal.N_GROUP_STRIP, abs=1e-4)
    assert float(lc.strip_index(0.40, 1.290)[0]) == pytest.approx(2.6192, abs=2e-3)
    assert float(lc.strip_index(0.40, 1.330)[0]) == pytest.approx(2.5672, abs=2e-3)
    widths = np.linspace(0.40, 1.0, 13)
    n_eff = [float(lc.strip_index(w)[0]) for w in widths]
    assert all(b > a for a, b in zip(n_eff, n_eff[1:])), "a wider strip has a higher index"
    with pytest.raises(ValueError):
        lc.strip_index(0.3)


# ---------------------------------------------------------------------------------------
# pin_attenuator
# ---------------------------------------------------------------------------------------
PIN_CASES = {
    "default": dict(),
    "short_narrow": dict(length=100.0, intrinsic_width=1.0, slab_taper=3.0, current=5.0),
    "long_wide": dict(length=500.0, intrinsic_width=2.4, slab_taper=6.0, current=10.0),
}


@pytest.mark.parametrize("case", sorted(PIN_CASES))
def test_pin_attenuator_interface(case):
    kw = PIN_CASES[case]
    cell = lc.pin_attenuator(**kw)
    length = kw.get("length", 200.0)
    assert sorted(cell.ports) == ["P0", "P1"] and sorted(cell.terminals) == ["ANODE", "CATHODE"]
    on_the_edge_and_pointing_in(cell)
    assert cell.ports["P0"].center[0] == pytest.approx(0.0, abs=1e-6) and cell.ports["P1"].center[0] == pytest.approx(length, abs=1e-6)
    # the ANODE (p) terminal is on the -y side, as on pn_phase_shifter
    assert cell.terminals["ANODE"].center()[1] < 0 < cell.terminals["CATHODE"].center()[1]
    # the slab is exactly the rib spec's 3.0 um at both port planes
    slab = np.asarray(cell.get_structures()[(2, 0)][0].vertices)
    for x in (0.0, length):
        ys = slab[np.isclose(slab[:, 0], x)][:, 1]
        assert ys.max() - ys.min() == pytest.approx(3.0)


@needs_klayout
@pytest.mark.parametrize("case", sorted(PIN_CASES))
def test_pin_attenuator_is_drc_clean(case):
    assert drc_clean(lc.pin_attenuator(**PIN_CASES[case])) == {}


@needs_klayout
def test_pin_attenuator_devrec_stops_at_the_ports_and_chains():
    """Two attenuators butted end to end are DEV.1 clean."""
    cell = lc.pin_attenuator(length=50.0)
    (dx0, _), (dx1, _) = cell.get_structures()[(100, 0)][0].bounds()
    assert dx0 == pytest.approx(0.0, abs=1e-6) and dx1 == pytest.approx(50.0, abs=1e-6)
    assert drc_clean(cell, [pf.Reference(lc.pin_attenuator(length=50.0), origin=(50.0, 0.0))]) == {}


def test_pin_attenuator_attenuation_follows_the_current_and_the_model_agrees():
    base = lc.pin_attenuator(length=500.0)
    t_base = float(np.abs(base.s_matrix(pf.C_0 / np.array([1.31]))[("P0@0", "P1@0")][0]) ** 2)
    last = 0.0
    for current in (1.0, 3.0, 10.0):
        cell = lc.pin_attenuator(length=500.0, current=current)
        atten = cell.properties["attenuation_db"]
        assert atten > last                                                      # monotonic in current
        last = atten
        t = float(np.abs(cell.s_matrix(pf.C_0 / np.array([1.31]))[("P0@0", "P1@0")][0]) ** 2)
        assert -10 * np.log10(t / t_base) == pytest.approx(atten, abs=0.01)     # the S matrix carries it
    assert 8 < lc.pin_attenuation_db(length=500.0, intrinsic_width=1.6, current=10.0) < 30   # a VOA, 10 mA in half a mm
    assert lc.pin_attenuation_db(length=500.0, intrinsic_width=1.6, current=0.0) == 0.0
    # the same current in a shorter device is a higher density, so the loss per mA rises as the device shrinks
    short = lc.pin_attenuation_db(length=100.0, intrinsic_width=1.6, current=5.0)
    long = lc.pin_attenuation_db(length=500.0, intrinsic_width=1.6, current=5.0)
    assert short > 0.5 * long


def test_pin_attenuator_refuses_a_current_beyond_the_range_of_the_absorption_fit():
    """The free-carrier fit ends at 1e20 cm^-3; 20 mA in 100 um x 1.0 um is 8.5e19 x 2 and raises, naming the limit."""
    with pytest.raises(ValueError, match="1e\\+20"):
        lc.pin_attenuation_db(length=100.0, intrinsic_width=1.0, current=20.0)
    with pytest.raises(ValueError, match="cm\\^-3"):
        lc.pin_attenuator(length=100.0, intrinsic_width=1.0, current=20.0)
    assert lc.pin_attenuation_db(length=100.0, intrinsic_width=1.0, current=5.0) > 0     # 4e19, inside the range


def test_pin_attenuator_names_follow_every_parameter_and_bad_input_raises():
    variants = [dict(), dict(length=100.0), dict(intrinsic_width=2.0), dict(slab_taper=5.0), dict(current=2.0)]
    names = [lc.pin_attenuator(**kw).name for kw in variants]
    assert len(set(names)) == len(names), names
    assert names[0] == "pin_l200p0"
    with pytest.raises(ValueError):
        lc.pin_attenuator(length=10.0)
    with pytest.raises(ValueError):
        lc.pin_attenuator(intrinsic_width=0.5)


# ---------------------------------------------------------------------------------------
# grating_antenna
# ---------------------------------------------------------------------------------------
ANT_CASES = {
    "default": dict(),
    "five_degrees": dict(period=lc.grating_antenna_period(5.0, width=0.8, duty=0.7), duty=0.7, length=60.0, width=0.8),
    "wide_period": dict(period=0.6, duty=0.6, length=200.0, width=0.5, taper_length=15.0),
}


@pytest.mark.parametrize("case", sorted(ANT_CASES))
def test_grating_antenna_interface_and_tooth_count(case):
    kw = ANT_CASES[case]
    cell = lc.grating_antenna(**kw)
    assert sorted(cell.ports) == ["P0"] and len(cell.terminals) == 0
    assert cell.ports["P0"].center == pytest.approx((0.0, 0.0)) and cell.ports["P0"].input_direction == 0.0
    n_teeth = len(cell.get_structures()[(3, 0)])
    assert n_teeth == int(kw.get("length", 100.0) // kw.get("period", 0.49))
    # every trench is as wide as period * (1 - duty) and as tall as the guide
    trench = [s.bounds() for s in cell.get_structures()[(3, 0)]][0]
    assert trench[1][0] - trench[0][0] == pytest.approx(kw.get("period", 0.49) * (1 - kw.get("duty", 0.75)), abs=1e-6)


@needs_klayout
@pytest.mark.parametrize("case", sorted(ANT_CASES))
def test_grating_antenna_is_drc_clean(case):
    assert drc_clean(lc.grating_antenna(**ANT_CASES[case])) == {}


@pytest.mark.parametrize("angle", [3.0, 10.0, 20.0])
def test_grating_antenna_period_solves_for_the_emission_angle(angle):
    period = lc.grating_antenna_period(angle)
    cell = lc.grating_antenna(period=period)
    assert cell.properties["emission_angle_deg"] == pytest.approx(angle, abs=0.3)      # a 1 nm period grid is 0.2 degrees


def test_grating_antenna_radiates_more_with_length_and_has_a_reflecting_port_model():
    fractions = [lc.grating_antenna(length=n).properties["radiated_fraction"] for n in (25.0, 50.0, 100.0, 200.0)]
    assert all(b > a for a, b in zip(fractions, fractions[1:])) and 0.0 < fractions[0] and fractions[-1] < 1.0
    # a longer antenna radiates from a longer aperture, so a narrower beam
    assert lc.grating_antenna(length=200.0).properties["divergence_deg"] < lc.grating_antenna(length=25.0).properties["divergence_deg"]
    s11 = np.abs(lc.grating_antenna().s_matrix(pf.C_0 / np.array([1.31]))[("P0@0", "P0@0")][0])
    assert 10 * np.log10(s11**2) == pytest.approx(-25.0, abs=0.1)


def test_grating_antenna_names_follow_every_parameter_and_deck_violations_raise():
    variants = [dict(), dict(period=0.5), dict(duty=0.7), dict(length=50.0), dict(width=0.8), dict(taper_length=12.0),
                dict(tip_length=8.0)]
    names = [lc.grating_antenna(**kw).name for kw in variants]
    assert len(set(names)) == len(names), names
    with pytest.raises(ValueError, match="GRAT"):
        lc.grating_antenna(period=0.49, duty=0.85)       # a 74 nm trench
    with pytest.raises(ValueError, match="GRAT"):
        lc.grating_antenna(period=0.49, duty=0.15)       # a 74 nm tooth
    with pytest.raises(ValueError):
        lc.grating_antenna(width=1.5)


# ---------------------------------------------------------------------------------------
# tw_electrode
# ---------------------------------------------------------------------------------------
TWE_CASES = {
    "default": dict(),
    "small": dict(length=100.0, signal_width=4.0, gap=4.0, ground_width=10.0, tap_pitch=5.0, pad_size=60.0, taper_length=60.0),
    "one_arm": dict(length=300.0, arms=1, ground_width=15.0, gap=8.0),
}


def _nets(cell):
    """Union-find over the conductive layers M1, VIA3 and M2 (M1-VIA3-M2) and VIA1 landing on M1, by shapely overlap."""
    import shapely
    from shapely.geometry import Polygon as SPoly, box

    def shape(s):
        if type(s).__name__ == "Rectangle" and abs(s.rotation) < 1e-9:
            (x0, y0), (x1, y1) = s.bounds()
            return box(x0, y0, x1, y1)
        return SPoly(np.asarray((s if type(s).__name__ == "Polygon" else s.to_polygon()).vertices))

    st = cell.get_structures()
    tech = cell.technology
    items = [(name, shape(s)) for name in ("M1", "VIA3", "M2", "VIA1") for s in st.get(tuple(tech.layers[name].layer), [])]
    parent = list(range(len(items)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    joins = {("M1", "VIA3"), ("VIA3", "M2"), ("M1", "VIA1")}
    tree = shapely.STRtree([g for _, g in items])
    for i, (n, g) in enumerate(items):
        for j in tree.query(g):
            m = items[j][0]
            if j > i and (n == m or (n, m) in joins or (m, n) in joins) and g.intersects(items[j][1]):
                parent[find(i)] = find(j)

    def net(layer, point):
        p = shapely.Point(point)
        for i, (n, g) in enumerate(items):
            if n == layer and g.contains(p):
                return find(i)
        return None

    return net


@pytest.mark.parametrize("case", sorted(TWE_CASES))
def test_tw_electrode_terminals_and_no_devrec(case):
    cell = lc.tw_electrode(**TWE_CASES[case])
    assert sorted(cell.terminals) == ["GAIN", "GAOUT", "GBIN", "GBOUT", "SIN", "SOUT"] and cell.ports == {}
    assert (100, 0) not in cell.get_structures()          # wiring: no DEVREC
    L, pm = TWE_CASES[case].get("length", 500.0), TWE_CASES[case].get("pad_size", 70.0) + 5.0
    T = TWE_CASES[case].get("taper_length", 100.0)
    sin, sout = cell.terminals["SIN"].center(), cell.terminals["SOUT"].center()
    assert sin[0] == pytest.approx(-T - pm / 2) and sout[0] == pytest.approx(L + T + pm / 2) and sin[1] == sout[1] == 0.0
    pitch = TWE_CASES[case].get("pad_pitch", 100.0)
    assert cell.terminals["GBIN"].center()[1] - cell.terminals["GAIN"].center()[1] == pytest.approx(2 * pitch)
    assert 30 < cell.properties["z0_ohm"] < 150 and cell.properties["n_microwave"] == pytest.approx(np.sqrt(6.35))


@needs_klayout
@pytest.mark.parametrize("case", sorted(TWE_CASES))
def test_tw_electrode_over_its_arms_is_drc_clean(case):
    assert drc_clean(lc.tw_electrode_with_arms(**TWE_CASES[case], name="tw_chk")) == {}


@needs_klayout
def test_tw_electrode_over_mirrored_arms_is_drc_clean():
    assert drc_clean(lc.tw_electrode_with_arms(length=100.0, push_pull=False, name="tw_chk_cm")) == {}


def _pad_and_arm_nets(top):
    net = _nets(top)
    ids = {name: net("M2", t.center()) for name, t in top.terminals.items() if name in ("SIN", "SOUT", "GAIN", "GBIN", "GAOUT", "GBOUT")}
    arms = {name: net("M1", t.center()) for name, t in top.terminals.items() if name[:2] in ("A_", "B_")}
    return ids, arms


def test_tw_electrode_makes_exactly_two_nets_and_the_arms_are_push_pull():
    """Default drive: the signal rail meets the CATHODE of arm A and the ANODE of arm B, the grounds the other two."""
    top = lc.tw_electrode_with_arms(length=300.0, name="tw_net")
    ids, arms = _pad_and_arm_nets(top)
    assert None not in ids.values() and None not in arms.values()
    assert len({ids[n] for n in ("GAIN", "GBIN", "GAOUT", "GBOUT")}) == 1, "all four ground pads are one net"
    assert ids["SIN"] == ids["SOUT"] != ids["GAIN"], "one signal net, apart from the ground"
    assert arms["A_CATHODE"] == arms["B_ANODE"] == ids["SIN"]
    assert arms["A_ANODE"] == arms["B_CATHODE"] == ids["GAIN"]


def test_tw_electrode_common_mode_drive_puts_both_cathodes_on_the_signal_rail():
    top = lc.tw_electrode_with_arms(length=300.0, push_pull=False, name="tw_net_cm")
    ids, arms = _pad_and_arm_nets(top)
    assert ids["SIN"] == ids["SOUT"] != ids["GAIN"]
    assert arms["A_CATHODE"] == arms["B_CATHODE"] == ids["SIN"]
    assert arms["A_ANODE"] == arms["B_ANODE"] == ids["GAIN"]


def test_tw_electrode_without_the_bridge_leaves_the_two_grounds_apart():
    top = lc.tw_electrode_with_arms(length=300.0, bridge=False, name="tw_nobridge")
    net = _nets(top)
    e = lc.tw_electrode(length=300.0, bridge=False)
    ids = {name: net("M2", t.center()) for name, t in e.terminals.items()}
    assert ids["GAIN"] == ids["GAOUT"] and ids["GBIN"] == ids["GBOUT"] and ids["GAIN"] != ids["GBIN"]


def test_tw_electrode_names_follow_every_parameter_and_bad_input_raises():
    variants = [dict(), dict(length=300.0), dict(signal_width=8.0), dict(gap=8.0), dict(ground_width=15.0), dict(arms=1),
                dict(contact_offset=1.2), dict(tap_pitch=20.0), dict(pad_size=60.0), dict(pad_pitch=120.0),
                dict(taper_length=80.0), dict(bridge=False)]
    names = [lc.tw_electrode(**kw).name for kw in variants]
    assert len(set(names)) == len(names), names
    for bad in (dict(arms=3), dict(gap=3.0), dict(length=10.0), dict(pad_pitch=75.0)):
        with pytest.raises(ValueError):
            lc.tw_electrode(**bad)


# ---------------------------------------------------------------------------------------
# contra_dc
# ---------------------------------------------------------------------------------------
CDC_CASES = {
    "default": dict(num_periods=950),
    "apodized": dict(num_periods=500, apodization=0.6, corrugation_width=0.08, gap=0.2),
    "channel_1271": dict(period=lc.contra_dc_period(1.271, width2=0.58), width2=0.58, num_periods=400, duty=0.33),
}


@pytest.mark.parametrize("case", sorted(CDC_CASES))
def test_contra_dc_interface(case):
    cell = lc.contra_dc(**CDC_CASES[case])
    assert sorted(cell.ports) == ["P0", "P1", "P2", "P3"] and cell.terminals == {}
    on_the_edge_and_pointing_in(cell)
    p = {n: cell.ports[n] for n in cell.ports}
    assert p["P0"].center[1] == p["P1"].center[1] == 0.0 and p["P2"].center[1] == p["P3"].center[1] < -6.0
    assert p["P0"].input_direction == p["P2"].input_direction == 0.0 and p["P1"].input_direction == p["P3"].input_direction == 180.0
    (dx0, _), (dx1, _) = cell.get_structures()[(100, 0)][0].bounds()
    assert dx0 == pytest.approx(0.0, abs=1e-6) and dx1 == pytest.approx(p["P1"].center[0], abs=1e-6)


@needs_klayout
@pytest.mark.parametrize("case", sorted(CDC_CASES))
def test_contra_dc_is_drc_clean(case):
    assert drc_clean(lc.contra_dc(**CDC_CASES[case])) == {}


@pytest.mark.parametrize("wavelength", [1.271, 1.291, 1.311, 1.331])
def test_contra_dc_drops_the_channel_it_was_designed_for(wavelength):
    """The four O-band CWDM channels: the period from `contra_dc_period` puts the drop peak on the channel."""
    period = lc.contra_dc_period(wavelength)
    # duty 0.34: at 1271 nm the period is 0.468 um and 0.32 of it is a 150 nm tooth, on the limit
    cell = lc.contra_dc(period=period, duty=0.34, num_periods=int(round(470.0 / period)))
    lam = np.linspace(1.25, 1.37, 6001)
    s = cell.s_matrix(pf.C_0 / lam)
    drop = np.abs(s[("P0@0", "P2@0")]) ** 2
    through = np.abs(s[("P0@0", "P1@0")]) ** 2
    assert abs(lam[np.argmax(drop)] - wavelength) < 5e-4            # 0.5 nm
    assert abs(cell.properties["bragg_wavelength_um"] - wavelength) < 5e-4
    assert 10 * np.log10(through[np.argmax(drop)]) < -20 and drop.max() > 0.95
    assert (through + drop).max() <= 1 + 1e-6                        # passive
    assert np.abs(s[("P0@0", "P3@0")]).max() < 1e-9                  # no co-directional coupling
    for a, b in (("P0", "P2"), ("P1", "P3"), ("P0", "P1"), ("P2", "P3")):
        assert np.abs(s[(a + "@0", b + "@0")] - s[(b + "@0", a + "@0")]).max() < 1e-9      # reciprocal


def test_contra_dc_periods_of_the_four_channels_are_increasing_and_second_order():
    periods = [lc.contra_dc_period(w) for w in (1.271, 1.291, 1.311, 1.331)]
    assert periods == sorted(periods) and all(0.46 < p < 0.52 for p in periods)
    assert lc.contra_dc_period(1.311, order=1) == pytest.approx(periods[2] / 2, abs=1e-4)


def test_contra_dc_apodization_lowers_the_side_lobes():
    lam = np.linspace(1.29, 1.33, 8001)

    def worst_sidelobe(cell):
        drop = np.abs(cell.s_matrix(pf.C_0 / lam)[("P0@0", "P2@0")]) ** 2
        half = lam[drop > 0.5][[0, -1]]
        out = (lam < half[0] - 1.5e-3) | (lam > half[1] + 1.5e-3)
        return 10 * np.log10(drop[out].max())

    uniform = lc.contra_dc(num_periods=500, corrugation_width=0.08)
    apod = lc.contra_dc(num_periods=500, corrugation_width=0.08, apodization=0.75)
    assert worst_sidelobe(apod) < worst_sidelobe(uniform) - 3.0


def test_contra_dc_first_order_and_tight_gap_raise_with_the_rule_named():
    with pytest.raises(ValueError, match="SI.W.1"):
        lc.contra_dc(period=lc.contra_dc_period(1.311, order=1), order=1, duty=0.5)
    with pytest.raises(ValueError, match="SI.S.1"):
        lc.contra_dc(gap=0.10)
    with pytest.raises(ValueError):
        lc.contra_dc(width2=0.42)                      # within 50 nm of guide 1: they would couple co-directionally
    with pytest.raises(ValueError):
        lc.contra_dc(apodization=0.8, corrugation_width=0.05)       # the shortest tooth would be 10 nm


def test_contra_dc_names_follow_every_parameter():
    variants = [dict(), dict(period=0.49), dict(num_periods=500), dict(width2=0.58), dict(width1=0.45), dict(gap=0.2),
                dict(corrugation_width=0.07), dict(duty=0.33), dict(order=3, period=0.737, duty=0.22), dict(apodization=0.5),
                dict(end_length=25.0), dict(port_offset=8.0)]
    names = [lc.contra_dc(**{"num_periods": 60, **kw}).name for kw in variants]
    assert len(set(names)) == len(names), names
