"""SiPhoPDK circuit-level cells: rings, delay lines, attenuators, antennas and electrodes.

They follow the conventions of `cells.py` (one `@pf.parametric_component` per cell,
recognition layers from `add_recognition`, a model that runs offline, clean against
`rules/siphopdk.drc`).

    from siphopdk import ladder_cells
    ring = ladder_cells.add_drop_ring(ring_radius=5.0, gap_through=0.16, gap_drop=0.16)

Cells
-----
add_drop_ring   a circular ring between two straight buses: four optical ports, an analytic
                `pf.RingModel` with one coupling coefficient per bus, and an optional heater.
add_drop_ring_corner
                the same ring between a row bus and a column bus that cross beside it, for a
                ring crossbar: the crossing itself is `siphopdk.crossing`, placed at the origin.
spiral_delay    a rectangular double spiral of exact centreline length, both ports side by side
                on one edge, with an analytic waveguide model.
pin_attenuator  a rib waveguide with a lateral p-i-n junction for carrier-injection attenuation:
                p and n implants either side of an intrinsic core, M1 contacts, two terminals, and
                a free-carrier-absorption model of attenuation against injected current.
grating_antenna a weak shallow-etch grating on a 1 um guide that radiates a few degrees off
                vertical: one optical port, a termination model, and the emission angle and
                radiated fraction stored on the cell.
tw_electrode    a travelling-wave GSG electrode on M2 for a push-pull pair of pn_phase_shifter arms
                (the arms lie in the two gaps), with M1 jumpers and VIA3 to the arms' contacts and
                a G-S-G pad set at each end; the ground pads are one net.
contra_dc       a contra-directional coupler: two corrugated guides of different width, four
                ports, a second-order grating (a first-order one at 1310 nm is under the 150 nm
                limit of SI.W.1 / SI.S.1) and a coupled-mode S matrix.
"""

import numpy as np
import photonforge as pf
import photonforge.typing as pft

from . import _calibration as cal
from ._recognition import add_recognition
from .cells import RIB, STRIP, WIDE, _GUIDE, _layer, _named, _profile, _spec, _tech, _wg_model

__all__ = ["add_drop_ring", "add_drop_ring_corner", "spiral_delay", "path_length", "ring_coupling", "ADD_DROP_PORTS", "MIN_RING_GAP",
           "pin_attenuator", "grating_antenna", "tw_electrode", "contra_dc",
           "strip_index", "contra_dc_period", "grating_antenna_period", "pin_attenuation_db"]

# The role of each optical port of `add_drop_ring`. The order is not a free choice: it is the
# order `pf.RingModel` expects for a double-bus ring (its S matrix couples P0-P2 and P1-P3 as
# the drop paths and P0-P1, P2-P3 as the through paths), so the layout and the model agree.
ADD_DROP_PORTS = {"in": "P0", "through": "P1", "drop": "P2", "add": "P3"}

# [DECK] SI.S.1: the narrowest legal space between two silicon shapes, in um. It is also the
# narrowest legal ring-to-bus gap, because that gap is a silicon space.
MIN_RING_GAP = 0.15

# [PHYS] oxide cladding index at the design wavelength (td Palik_NoLoss SiO2, 1.31 um).
_N_CLAD = 1.447

# Where the heater leads leave the ring: the same 135 and 45 degrees as the modulator and the
# all-pass ring, so the two rings are interchangeable when routing heater metal.
_LEAD_DIRECTION, _LEAD_SPREAD = 90.0, 90.0


# =====================================================================
# Coupling model
# =====================================================================
def _field_decay(n_eff):
    """Decay constant of the guided field in the oxide, 1/um: k0 * sqrt(n_eff^2 - n_clad^2)."""
    return 2 * np.pi / cal.WAVELENGTH * np.sqrt(n_eff**2 - _N_CLAD**2)


def _effective_length(radius, coupling_length, gamma):
    """Length over which a straight bus and a ring of this radius really couple, um.

    The gap grows as the ring curves away, g(x) = g0 + x^2 / (2R), and the field coupling
    falls as exp(-gamma g), so the coupling integrates to a Gaussian of width sqrt(R / gamma):
    the integral is sqrt(2 pi R / gamma). Any straight run adds to it one for one.
    """
    return coupling_length + np.sqrt(2 * np.pi * radius / gamma)


def _anchor():
    """The one calibrated point: the ring that `cal.KAPPA1` was measured on.

    The calibration is the 12 um ring with a 180 nm gap and a 4 um straight run, coupling
    0.387. Everything else is scaled from it with the physics above, so the constants stay
    the PDK's own and nothing is fitted to a foundry.
    """
    gamma = _field_decay(cal.N_EFF_RIB)
    k_ref = np.arcsin(cal.KAPPA1)
    return k_ref / (np.exp(-gamma * 0.18) * _effective_length(cal.R_CH0, 4.0, gamma))


def ring_coupling(radius, gap, *, coupling_length=0.0, port_spec=STRIP):
    """Field coupling coefficient between a ring and a bus, from the geometry. [SYNTH]

    `radius` is the ring centreline radius and `gap` the edge-to-edge space, both in um. This
    is the number `pf.RingModel` takes as `kappa`. It is anchored to the calibrated ring in
    `_calibration.py` and moves with the gap as exp(-gamma g) and with the radius through the
    coupling length, so it is a trend with the right shape, not a prediction.
    """
    n_eff = _GUIDE[port_spec][0]
    gamma = _field_decay(n_eff)
    k = _anchor() * np.exp(-gamma * gap) * _effective_length(radius, coupling_length, gamma)
    return float(np.sin(min(k, 1.4)))  # a coupler cannot move more than all of the light


# =====================================================================
# Heater, shared by the ring cells
# =====================================================================
def _tip_radius_clear_of_buses(ring_radius, strip_width, y_thru, y_drop):
    """How far out along its lead the landing pad sits, as a function of the lead direction.

    R + 6 um is the modulator's, but on a small ring that puts the 5.2 um pad over the bus the
    lead has just crossed (at R = 5 the pad covers y = 5.2 to 10.4 um and bus 2 is at 5.6 um).
    The pad is pushed out until its nearer edge clears the far bus by 0.5 um, so a router has
    a clean landing; for rings of about 8 um and up the answer is still R + 6.
    """
    clear = 2.6 + 0.5 + strip_width / 2

    def tip_radius(u):
        r = ring_radius + 6.0
        if u[1] > 0.05:
            r = max(r, (y_drop + clear) / u[1])
        elif u[1] < -0.05:
            r = max(r, (-y_thru + clear) / -u[1])
        return r

    return tip_radius


def _draw_heater(c, tech, *, center, ring_radius, strip_width, heater_width, lead_direction,
                 lead_spread, tip_radius):
    """A TiN arc over the ring with a radial lead and an M1 landing pad at each end.

    The arc covers the whole ring except the wedge `lead_spread` wide centred on
    `lead_direction`; HTRP is the lead at the larger angle and HTRN the one at the smaller.
    Adds the HTR, VIA2 and M1 shapes and the two terminals to `c`.
    """
    cx, cy = center
    a_plus = lead_direction + lead_spread / 2
    a_minus = lead_direction - lead_spread / 2
    # Each lead starts on the ring centreline, in the middle of the arc's band. Starting it
    # nearer the ring's centre leaves a notch between the lead and the arc, and on a 3 um
    # ring that notch is under the 1 um HTR.S.1 space (the modulator's leads start nearer
    # the centre, and its 12 um ring has room for it). Starting it exactly at the inner
    # edge makes the lead's corner graze the arc's polygonal inner circle, which leaves a
    # 1 nm sliver at some radii (5, 6, 7 and 15 um). Starting 0.3 um or more into the band
    # was clean at every radius tried from 3 to 15 um. Two leads closer than the spacing
    # rule would merge or fail, and that shows up as a DRC marker far from the cause, so
    # check it here: the chord between the lead centrelines at their inner ends, less the
    # lead width, is the space where they are closest.
    inner = ring_radius
    if ring_radius - heater_width / 2 <= 0:
        raise ValueError(f"heater_width {heater_width} um is too wide for a {ring_radius} um ring")
    chord = 2 * inner * np.sin(np.radians(lead_spread) / 2) - heater_width
    if chord < 1.1:
        raise ValueError(
            f"the heater leads are {chord:.2f} um apart at their inner ends, under the "
            "1 um HTR.S.1 space; open lead_spread or use a bigger ring"
        )

    def unit(angle):
        return np.array([np.cos(np.radians(angle)), np.sin(np.radians(angle))])

    origin = np.array([cx, cy], float)
    c.add("HTR", pf.Circle(radius=ring_radius + heater_width / 2,
                           inner_radius=ring_radius - heater_width / 2, center=(cx, cy),
                           sector=(a_plus, a_minus + 360)))
    for angle in (a_plus, a_minus):
        u = unit(angle)
        n = np.array([-u[1], u[0]]) * (heater_width / 2)
        lo, hi = origin + u * inner, origin + u * (tip_radius(u) + 1.2)
        c.add("HTR", pf.Polygon([tuple(float(v) for v in lo - n), tuple(float(v) for v in hi - n),
                                 tuple(float(v) for v in hi + n), tuple(float(v) for v in lo + n)]))
    # HTRP and HTRN, not H+ and H-: the name becomes a net name in a SPICE netlist and `+` is an
    # operator there (see micro_ring_modulator). Pads are snapped to the grid and wider than
    # the 4 um trace that lands on them, for the same reason as the modulator's.
    leads = {}
    for name, angle in zip(("HTRP", "HTRN"), (a_plus, a_minus)):
        u = unit(angle)
        tip = tuple(round(float(v), 3) for v in origin + u * tip_radius(u))
        c.add("VIA2", pf.Rectangle(center=tip, size=(0.4, 0.4)))
        c.add("M1", pf.Rectangle(center=tip, size=(5.2, 5.2)))
        leads[name] = pf.Terminal(_layer(tech, "M1"), pf.Rectangle(center=tip, size=(5.2, 5.2)))
    c.add_terminal(leads)


# =====================================================================
# Add-drop ring
# =====================================================================
@pf.parametric_component(name_prefix="adr")
def add_drop_ring(*, ring_radius: pft.PositiveDimension = 5.0,
                  gap_through: pft.PositiveDimension = 0.18,
                  gap_drop: pft.PositiveDimension = 0.18,
                  heater: bool = True,
                  heater_width: pft.PositiveDimension = 1.2,
                  lead_direction: pft.Angle = _LEAD_DIRECTION,
                  lead_spread: pft.Angle = _LEAD_SPREAD,
                  bus_extension: pft.PositiveDimension = 6.0,
                  temperature: pft.Temperature = 293.0,
                  port_spec: str = STRIP, technology=None):
    """Add-drop microring: a circular ring between two straight buses, with an optional heater.

        P2 drop  <----------- bus 2, gap_drop ----------<  P3 add
                         .--------.
                        (   ring   )   heater arc and its two leads (HTRP, HTRN)
                         '--------'
        P0 in    >----------- bus 1, gap_through -------->  P1 through

    The ring is centred on the origin and its centreline radius is `ring_radius`. Bus 1 is the
    lower one and carries `in` and `through`; bus 2 is the upper one and carries `drop` and
    `add`. Light entering `in` runs counter-clockwise round the ring and leaves at `drop`.
    The gaps are edge to edge, and they are silicon spaces, so the deck's SI.S.1 sets a floor
    of 150 nm: a smaller gap raises here instead of failing DRC later.

    Port angles follow the PDK: the angle is the direction light travels into the cell, so
    every port sits on the outside of the cell and looks out along the bus. `bus_extension`
    is how far each bus runs past the ring's own width.

    The ring is a 400 nm strip, not the rib of the modulator. A filter has no junction, so
    it needs no slab; a strip bends tighter, which is what a 3 um ring asks for; and the
    grating couplers are strip too, so nothing needs a transition. Strip loss is in
    `_calibration.py`, and bend radiation is not modelled, so treat radii below about 5 um as
    a layout test and not as a loss prediction.

    The heater is an arc of TiN over the ring, open where the two leads leave. Leads cross
    bus 2 on the way out when `lead_direction` is 90 (up). They are over the waveguide, not
    in it, and that is the same trade the modulator makes; set `lead_direction` to 180 or 0
    to leave between the buses instead, when the ring is big enough to allow it.

    `temperature` is a cell parameter for the same reason `voltage` is on the modulator: the
    S matrix is cached per component, so changing it on the model afterwards returns the old
    answer.
    """
    if port_spec not in (STRIP, WIDE):
        raise ValueError(
            f"add_drop_ring is drawn on a single silicon layer, so port_spec must be "
            f"{STRIP!r} or {WIDE!r}, not {port_spec!r}: a rib needs a slab drawn round the ring "
            "and a nitride guide has its own rules"
        )
    for name, gap in (("gap_through", gap_through), ("gap_drop", gap_drop)):
        if gap < MIN_RING_GAP - 1e-9:
            raise ValueError(
                f"{name} = {gap * 1e3:.0f} nm is below the 150 nm silicon space of SI.S.1. "
                "The ring-to-bus gap is a silicon space, so the deck forbids it."
            )

    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    # ring, centred on the origin
    c.add("SI", pf.Circle(radius=ring_radius + w / 2, inner_radius=ring_radius - w / 2))

    # two straight buses tangent to the ring, one below and one above. The centreline of a bus
    # is ring_radius + w (half a ring width plus half a bus width) plus the gap from the ring.
    half = ring_radius + w + bus_extension
    y_thru = -(ring_radius + w + gap_through)
    y_drop = ring_radius + w + gap_drop
    for y in (y_thru, y_drop):
        c.add("SI", pf.Path((-half, y), w).segment((2 * half, 0), relative=True))

    # the heater: a broken arc over the ring, one radial lead at each end. The landing pad sits
    # along its lead far enough out to clear the far bus (see _tip_radius_clear_of_buses).
    if heater:
        _draw_heater(c, tech, center=(0.0, 0.0), ring_radius=ring_radius, strip_width=w,
                     heater_width=heater_width, lead_direction=lead_direction,
                     lead_spread=lead_spread,
                     tip_radius=_tip_radius_clear_of_buses(ring_radius, w, y_thru, y_drop))

    # Four ports, in the order pf.RingModel wants (ADD_DROP_PORTS). The port angle is the
    # direction of travel into the cell, so the left ones look right and the right ones left.
    c.add_port(pf.Port((-half, y_thru), 0, spec), "P0")
    c.add_port(pf.Port((half, y_thru), 180, spec), "P1")
    c.add_port(pf.Port((-half, y_drop), 0, spec), "P2")
    c.add_port(pf.Port((half, y_drop), 180, spec), "P3")

    # A circle, so the round trip is exactly 2 pi R and the index belongs to the cross-section.
    n_eff, n_group, loss = _GUIDE[port_spec]
    length = 2 * np.pi * ring_radius
    kappa1 = ring_coupling(ring_radius, gap_through, port_spec=port_spec)
    kappa2 = ring_coupling(ring_radius, gap_drop, port_spec=port_spec)
    model = pf.RingModel(kappa1=kappa1, kappa2=kappa2, n_eff=n_eff, n_group=n_group,
                         length=length, propagation_loss=loss,
                         reference_frequency=cal.REF_FREQUENCY,
                         dn_dT=cal.DN_DT, temperature=temperature)
    model.time_stepper = pf.RingTimeStepper(
        kappa1=kappa1, kappa2=kappa2, n_eff=n_eff, n_group=n_group, length=length,
        propagation_loss=loss, dn_dT=cal.DN_DT, temperature=temperature)
    c.add_model(model, "ring")
    c.activate_model("ring")

    # The name is the cell's identity wherever a name is the key (the viewer palette, PDA, a GDS
    # written and read back), so it carries every parameter that changes the cell. The first four
    # are always there, which keeps the names that already existed; the rest are added only when
    # they differ from the default, so those names did not change either.
    extras = {}
    for key, value, default in (("be", bus_extension, 6.0), ("hw", heater_width, 1.2),
                                ("ld", lead_direction, _LEAD_DIRECTION), ("ls", lead_spread, _LEAD_SPREAD),
                                ("t", temperature, 293.0), ("ps", port_spec, STRIP)):
        if value != default:
            extras[key] = value
    _named(c, "adr", r=ring_radius, gt=gap_through, gd=gap_drop, h=int(heater), **extras)

    # DEVREC ends exactly at the port plane, like the PDK's straight and termination. The default
    # box comes from c.bounds(), which ran 0.08 um past the ports (the bus ends are +-half but the
    # bounds said +-(half + 0.08)), so two of these butted together, or one butted to a straight,
    # overlapped by 0.16 um and DEV.1 reported two devices on the same spot. The box is built from
    # the drawn shapes instead, with x held to the ports.
    boxes = np.array([s.bounds() for shapes in c.get_structures().values() for s in shapes])
    y0, y1 = float(boxes[:, 0, 1].min()), float(boxes[:, 1, 1].max())
    add_recognition(c, icon="mrr", devrec=False)
    c.add("DEVREC", pf.Rectangle((-half, y0), (half, y1)))
    return c


# =====================================================================
# Add-drop ring at the corner of a crossing: the crosspoint of a ring crossbar
# =====================================================================
@pf.parametric_component(name_prefix="adrc")
def add_drop_ring_corner(*, ring_radius: pft.PositiveDimension = 10.0,
                         gap_row: pft.PositiveDimension = 0.18,
                         gap_column: pft.PositiveDimension = 0.18,
                         heater: bool = True,
                         heater_width: pft.PositiveDimension = 1.2,
                         lead_direction: pft.Angle = 45.0,
                         lead_spread: pft.Angle = 90.0,
                         arm_length: pft.PositiveDimension = 6.0,
                         extension: pft.PositiveDimension = 3.0,
                         temperature: pft.Temperature = 293.0,
                         port_spec: str = STRIP, technology=None):
    """Add-drop microring between a row bus and a column bus that cross next to it.

        P3 column top
         |
         |     .--------.
         |    (   ring   )   heater arc, leads to HTRP (up) and HTRN (right)
         |     '--------'
         |  ______________________
        (x)-/  P0 row start --->  P1 row end
         |
        P2 column start          (x) = the centre of a waveguide crossing, at the origin

    The origin is the centre of the crossing this cell is meant to sit beside, and the cell
    does not draw it: place `siphopdk.crossing(arm_length=arm_length)` at the origin and butt
    the cell's P0 and P2 to the crossing's right and top arms. The row bus then runs from the
    crossing's right arm, under the ring, to P1; the column bus runs from its top arm, past
    the ring's left side, to P3. The ring is in the upper-right quadrant, tangent to both.

    Ports follow `ADD_DROP_PORTS`, with the row as bus 1 and the column as bus 2, so the
    `pf.RingModel` ordering and the layout agree: P0 row start, P1 row end, P2 column start,
    P3 column top. Light entering P0 (moving right) runs counter-clockwise round the ring and
    leaves down the column at P2; light entering P1 leaves up the column at P3. A signal that
    comes down the column from P3 is dropped to the row at P1. Port angles follow the PDK.

    The two leads of the heater exit at 0 and 90 degrees, into the empty quadrant, so the
    landing pads touch neither bus and no heater metal sits over a waveguide. A ring
    under 6 um raises, because it would sit within 150 nm of the crossing's tapered arms (DRC
    flagged 8 SI.S.1 at 5 um); from 6 um up the cell plus a crossing is DRC clean. Keep it at
    8 um or more if the gap to the arms should stay well above the coupling gap. Gaps below
    150 nm raise, as in `add_drop_ring`.
    """
    if port_spec not in (STRIP, WIDE):
        raise ValueError(f"add_drop_ring_corner is drawn on one silicon layer: port_spec must be "
                         f"{STRIP!r} or {WIDE!r}, not {port_spec!r}")
    if ring_radius < 6.0:
        raise ValueError(f"ring_radius {ring_radius} um is too small for add_drop_ring_corner: below "
                         "6 um the ring is within 150 nm of the crossing's tapered arms (SI.S.1)")
    for name, gap in (("gap_row", gap_row), ("gap_column", gap_column)):
        if gap < MIN_RING_GAP - 1e-9:
            raise ValueError(f"{name} = {gap * 1e3:.0f} nm is below the 150 nm silicon space of "
                             "SI.S.1. The ring-to-bus gap is a silicon space, so the deck forbids it.")

    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    # ring centre: tangent to the row bus (y = 0) below it and the column bus (x = 0) beside it
    cx = ring_radius + w + gap_column
    cy = ring_radius + w + gap_row
    c.add("SI", pf.Circle(radius=ring_radius + w / 2, inner_radius=ring_radius - w / 2, center=(cx, cy)))
    x_end = cx + ring_radius + extension
    y_end = cy + ring_radius + extension
    c.add("SI", pf.Path((arm_length, 0), w).segment((x_end - arm_length, 0), relative=True),
          pf.Path((0, arm_length), w).segment((0, y_end - arm_length), relative=True))

    if heater:
        _draw_heater(c, tech, center=(cx, cy), ring_radius=ring_radius, strip_width=w,
                     heater_width=heater_width, lead_direction=lead_direction,
                     lead_spread=lead_spread, tip_radius=lambda u: ring_radius + 6.0)

    c.add_port(pf.Port((arm_length, 0), 0, spec), "P0")      # row start, looks along +x
    c.add_port(pf.Port((x_end, 0), 180, spec), "P1")         # row end
    c.add_port(pf.Port((0, arm_length), 90, spec), "P2")     # column start, looks along +y
    c.add_port(pf.Port((0, y_end), -90, spec), "P3")         # column top

    n_eff, n_group, loss = _GUIDE[port_spec]
    length = 2 * np.pi * ring_radius
    kappa1 = ring_coupling(ring_radius, gap_row, port_spec=port_spec)
    kappa2 = ring_coupling(ring_radius, gap_column, port_spec=port_spec)
    model = pf.RingModel(kappa1=kappa1, kappa2=kappa2, n_eff=n_eff, n_group=n_group,
                         length=length, propagation_loss=loss,
                         reference_frequency=cal.REF_FREQUENCY,
                         dn_dT=cal.DN_DT, temperature=temperature)
    model.time_stepper = pf.RingTimeStepper(
        kappa1=kappa1, kappa2=kappa2, n_eff=n_eff, n_group=n_group, length=length,
        propagation_loss=loss, dn_dT=cal.DN_DT, temperature=temperature)
    c.add_model(model, "ring")
    c.activate_model("ring")

    extras = {}
    for key, value, default in (("hw", heater_width, 1.2), ("ld", lead_direction, 45.0),
                                ("ls", lead_spread, 90.0), ("al", arm_length, 6.0),
                                ("ex", extension, 3.0), ("t", temperature, 293.0), ("ps", port_spec, STRIP)):
        if value != default:
            extras[key] = value
    _named(c, "adrc", r=ring_radius, gr=gap_row, gc=gap_column, h=int(heater), **extras)

    # DEVREC is the device boundary, normally the cell's bounding box. Here the box would run
    # into the crossing this cell is butted to (its own boundary is the square of half-width
    # `arm_length`) and DEV.1 would read two overlapping devices, so the boundary is the box with
    # that quadrant cut away: an L that touches the crossing's edges and does not overlap them.
    # The extent comes from the drawn shapes, not c.bounds(), which includes the 2.4 um width of
    # the ports and would reach back into the crossing.
    boxes = np.array([s.bounds() for shapes in c.get_structures().values() for s in shapes])
    (x0, y0), (x1, y1) = boxes[:, 0].min(axis=0), boxes[:, 1].max(axis=0)
    a = arm_length
    c.add("DEVREC", pf.Polygon([(a, y0), (x1, y0), (x1, y1), (x0, y1), (x0, a), (a, a)]))
    return add_recognition(c, icon="mrr", devrec=False)


# =====================================================================
# Spiral delay line
# =====================================================================
# spiral_delay: a rectangular double spiral ("paperclip") delay line for SiPhoPDK.
# 
# The spiral is written with PhotonForge paths only (segments and turns) and follows the SiEPIC `spiral_paperclip` topology that the openEBL MZI
# reference uses: two interleaved arms wind inwards together, meet in an S-bend at the centre,
# and the path comes back out along the second arm, so BOTH ports sit side by side on the same
# (left) edge, `separation` apart, both facing left.
# 
#     P0 -->--------------------------------------+   outer arm, winds in (clockwise)
#     P1 -->------------------------------+       |
#          |  +---------------------------+ |     |
#          |  |  +---------------------+  | |     |
#          |  |  |        S            |  | |     |   the centre S-bend joins the two arms
#          |  |  +----------------------+ | |     |
#          |  +------------------------------+    |
#          +--------------------------------------+
# 
# The path length is exact by construction: every lane length is linear in the cell width, so
# the width is solved so that the true centreline length equals `length`. The tests check the
# result against `photonforge.route_length`.
# 
# Geometry (all centreline coordinates, y down from the top lane, x right from the ports):
# 
#     ring j = 0 .. 2n-1, arm A owns even rings and arm B odd rings
#     top lane      y_t(j) = -p j              (p = separation)
#     bottom lane   y_b(j) = -H + p j          (H = 2 R + p (4 n - 3), so the centre gap is 2 R)
#     right link    x_r(j) = W - p j           (j = 0 .. 2n-2)
#     left link     x_l(j) = lead + p j        (j = 0 .. 2n-3)
# 
# Every corner is a 90 degree `Path.turn` of radius R (the PDK bend: Euler fraction 0.2), whose
# end point is exactly R away on each axis, so a corner needs R of straight on each side.

_EPS = 1e-9


# ---------------------------------------------------------------------------------------
# polyline -> fillet operations
# ---------------------------------------------------------------------------------------
def _unit(a, b):
    d = np.asarray(b, float) - np.asarray(a, float)
    n = np.hypot(*d)
    return d / n, n


def _fillet_ops(vertices, radius):
    """Turn a polyline with right-angle corners into [("seg", dx, dy) | ("turn", +-90)].

    A corner of radius R consumes R of each neighbouring segment. Returns the operations and
    the list of straight lengths, so a caller can reject an infeasible spiral.
    """
    v = [np.asarray(p, float) for p in vertices]
    ops, straights = [], []
    pos = v[0].copy()
    for k in range(1, len(v)):
        d_in, _ = _unit(v[k - 1], v[k])
        if k < len(v) - 1:
            d_out, _ = _unit(v[k], v[k + 1])
            cross = d_in[0] * d_out[1] - d_in[1] * d_out[0]
            if abs(cross) < 0.5:
                raise ValueError("spiral_delay: only right-angle corners are supported")
            tangent = v[k] - d_in * radius        # where the arc starts
            gap = tangent - pos
            straight = float(np.dot(gap, d_in))   # signed length along the heading
            straights.append(straight)
            if straight > _EPS:
                ops.append(("seg", float(gap[0]), float(gap[1])))
            ops.append(("turn", 90.0 if cross > 0 else -90.0))
            pos = v[k] + d_out * radius           # where the arc ends
        else:
            gap = v[k] - pos
            straight = float(np.dot(gap, d_in))
            straights.append(straight)
            if straight > _EPS:
                ops.append(("seg", float(gap[0]), float(gap[1])))
            pos = v[k]
    return ops, straights


def _spiral_vertices(turns, p, radius, width, lead):
    """Sharp-corner centreline of the two-arm spiral (see the module docstring)."""
    n = turns
    h = 2 * radius + p * (4 * n - 3)
    yt = lambda j: -p * j
    yb = lambda j: -h + p * j
    xr = lambda j: width - p * j
    xl = lambda j: lead + p * j
    x_c = 0.5 * (xr(2 * n - 2) + xl(2 * n - 3)) if n > 1 else 0.5 * (xr(0) + lead)

    pts = [(0.0, yt(0))]
    # arm A: inwards, clockwise
    for k in range(n):
        j = 2 * k
        pts += [(xr(j), yt(j)), (xr(j), yb(j))]
        if j < 2 * n - 2:
            pts += [(xl(j), yb(j)), (xl(j), yt(j + 2))]
        else:
            pts += [(x_c, yb(j)), (x_c, yt(2 * n - 1))]
    # arm B: outwards, counter-clockwise (ring 2n-1 is only the centre top lane)
    for k in range(n - 1, 0, -1):
        j = 2 * k - 1
        pts += [(xl(j), yt(j + 2)), (xl(j), yb(j)), (xr(j), yb(j)), (xr(j), yt(j))]
    pts += [(0.0, yt(1))]
    # drop the duplicated corner that appears when a link has zero straight length
    out = [pts[0]]
    for q in pts[1:]:
        if np.hypot(q[0] - out[-1][0], q[1] - out[-1][1]) > _EPS:
            out.append(q)
    return out, h


def _bend_probe(spec_width, radius, euler_fraction):
    """Length of one 90 degree turn (the corner footprint is exactly `radius`)."""
    path = pf.Path((0, 0), spec_width).turn(-90, radius=radius, euler_fraction=euler_fraction)
    end = path.spine()[-1]
    assert abs(abs(end[0]) - radius) < 1e-6 and abs(abs(end[1]) - radius) < 1e-6, end
    return float(path.length())


def analytic_length(vertices, radius, bend_len):
    """Sum of straight pieces plus one bend per corner (what the path should measure)."""
    _, straights = _fillet_ops(vertices, radius)
    return float(sum(straights) + bend_len * (len(vertices) - 2))


def _solve(length, turns, p, radius, lead, bend_len):
    """(width, height, vertices, dL/dW) that make the centreline `length`, or None if it does not fit."""
    def total(w):
        vs, _ = _spiral_vertices(turns, p, radius, w, lead)
        return analytic_length(vs, radius, bend_len)

    # every lane is linear in W, so two probes fix the line L = a + b W exactly
    w0, w1 = 100.0 * radius + 20.0 * p * turns, 200.0 * radius + 20.0 * p * turns
    l0, l1 = total(w0), total(w1)
    slope = (l1 - l0) / (w1 - w0)
    width = w0 + (length - l0) / slope
    vs, h = _spiral_vertices(turns, p, radius, width, lead)
    _, straights = _fillet_ops(vs, radius)
    if width <= 0 or min(straights) < -1e-6:
        return None
    return width, h, vs, slope


@pf.parametric_component(name_prefix="spiral")
def spiral_delay(*, length: pft.PositiveDimension = 4400.0,
                 separation: pft.PositiveDimension = 4.0,
                 radius: pft.Dimension = 0.0,
                 aspect: pft.PositiveDimension = 8.0,
                 turns: int = 0,
                 lead: pft.PositiveDimension = 1.0,
                 port_spec: str = STRIP,
                 euler_fraction: float = 0.2,
                 technology=None):
    """Rectangular double-spiral delay line with both ports on the left edge.

    Args:
        length: true centreline length in um (met to better than a nanometre).
        separation: lane pitch p in um. Two lanes of the same arm are 2 p apart.
        radius: corner radius; 0 means the port spec's `default_radius`.
        aspect: target width / height of the finished spiral; ignored when `turns` is set.
        turns: turns per arm; 0 picks the count whose shape is closest to `aspect`.
        lead: x of the first left-hand link, so no silicon sticks out past the port plane.
        port_spec: guide the spiral is drawn in (SI_STRIP_O, or SIN_STRIP_O for nitride).

    The cell carries `properties["path_length_um"]` (true length) and an analytic model.
    Ports P0 (outer arm, top) and P1 (inner arm, `separation` below it) both face left.
    """
    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    r = radius if radius > 0 else spec.default_radius
    if r <= 0:
        raise ValueError(f"spiral_delay: {port_spec!r} has no default_radius, pass radius=")
    w_line = spec.path_profile_for(_layer(tech, "SIN" if port_spec.startswith("SIN") else "SI"), tech)[0]
    bend_len = _bend_probe(w_line, r, euler_fraction)

    # choose the turn count: width falls and height grows with n, so walk up until the shape
    # passes the requested aspect ratio and keep whichever neighbour is closer to it
    best = None
    for n in ([turns] if turns > 0 else range(1, 400)):
        sol = _solve(length, n, separation, r, lead, bend_len)
        if sol is None:
            continue
        width, h = sol[0], sol[1]
        score = abs(np.log(width / h / aspect))
        if best is None or score < best[0]:
            best = (score, n, sol)
        if turns <= 0 and width / h < aspect:
            break
    if best is None:
        raise ValueError(
            f"spiral_delay: {length} um does not fit with radius {r} and separation "
            f"{separation} (too short for the bends, or no feasible turn count)")
    _, n, (width, height, vertices, slope) = best

    # build, measure the real path, and correct the width once for the bend-length rounding
    for attempt in range(2):
        ops, _ = _fillet_ops(vertices, r)
        c = pf.Component(technology=tech)
        path_len = None
        for layer, path in spec.get_paths((0, 0)):
            for op in ops:
                if op[0] == "seg":
                    path.segment((op[1], op[2]), relative=True)
                else:
                    path.turn(op[1], radius=r, euler_fraction=euler_fraction)
            c.add(layer, path)
            if path_len is None:
                path_len = float(path.length())
        if attempt == 0 and abs(path_len - length) > 1e-5:
            width += (length - path_len) / slope
            vertices, height = _spiral_vertices(n, separation, r, width, lead)
        else:
            break

    c.add_port(pf.Port((0.0, 0.0), 0, spec), "P0")
    c.add_port(pf.Port((0.0, -separation), 0, spec), "P1")

    n_bends = len(vertices) - 2
    c.add_model(_wg_model(path_len, port_spec=port_spec,
                                 extra_loss=cal.LOSS_BEND_90 * n_bends), "analytic")
    c.activate_model("analytic")
    c.properties["path_length_um"] = path_len
    c.properties["turns_per_arm"] = n
    # The name is the cell's identity wherever a name is the key, so it carries every parameter that
    # changes the geometry. The first four are always there; the rest only when not the default.
    extras = {}
    for key, value, default in (("ef", euler_fraction, 0.2), ("ld", lead, 1.0), ("ps", port_spec, STRIP)):
        if value != default:
            extras[key] = value
    _named(c, "spiral", L=length, p=separation, r=r, n=n, **extras)
    # DEVREC hugs the cell like every PDK cell, except that the ports sit on x = 0 and a
    # bounding box taken with them is 1 nm wider, which makes a butt-connected transition
    # overlap it by 1 nm and trips DEV.1. So the box starts exactly at the port plane.
    (_, y0), (x1, y1) = c.bounds()
    add_recognition(c, icon="spiral", devrec=False)
    c.add("DEVREC", pf.Rectangle((0.0, y0), (x1, y1)))
    return c


def path_length(component):
    """True centreline length stored by `spiral_delay` (falls back to `pf.route_length`)."""
    try:
        return float(component.properties["path_length_um"])
    except KeyError:
        return float(pf.route_length(component))


# =====================================================================
# pin_attenuator, grating_antenna, tw_electrode, contra_dc
# =====================================================================
# The same rules as the cells above: the name
# carries every parameter, DEVREC ends at the port plane (a metal-only cell has none), every
# layer tuple comes through the cell's own technology, and a model runs offline. Constants are
# tagged [SOLVED] (a local mode solve on this stack), [PHYS] (published physics), [DECK] (a rule
# of rules/siphopdk.drc) or [SYNTH] (invented, plausible, easy to replace).

_GDS_MARGIN = 0.003        # [DECK] a feature drawn this much over a 150 nm rule still passes after the 1 nm grid snaps it
_RULE_SI = 0.150           # [DECK] SI.W.1 and SI.S.1, um
_RULE_GRAT = 0.100         # [DECK] GRAT.W.1 and GRAT.S.1, um

# ---------------------------------------------------------------------------------------
# Effective index of a silicon strip against its width
# ---------------------------------------------------------------------------------------
# [SOLVED] The local (free) mode solver on this PDK's own stack, 1310 nm, the same call as
# scripts/calibrate.py. The 400 nm row is the calibration table's SI_STRIP_O (2.5932, 4.2949) to four
# digits, so everything below is on the PDK's footing and not a separate guess. Widths other than
# these are interpolated linearly. The contra-DC needs two widths, and the table in
# _calibration.py has one strip.
_STRIP_WIDTHS = np.array([0.40, 0.46, 0.52, 0.58, 0.70, 1.00])
_STRIP_NEFF = np.array([2.5932, 2.6844, 2.7470, 2.7919, 2.8505, 2.9176])
_STRIP_NGROUP = np.array([4.2949, 4.1859, 4.1090, 4.0534, 3.9805, 3.8970])
# [SOLVED] Index lost by a strip when the silicon is 150 nm thick instead of 220 nm (the floor of a
# 70 nm shallow-etch trench): 1.0 um 2.9176 -> 2.6249, 0.7 um 2.8505 -> 2.5499, so 0.293 to 0.301.
_ETCH_DN = 0.295


def strip_index(width, wavelength=cal.WAVELENGTH):
    """(n_eff, n_group) of a 220 nm silicon strip of this width (um, 0.40 to 1.00) at `wavelength` (um).

    The width dependence is the [SOLVED] table above. The wavelength dependence is first order, n_eff
    = n0 + (n0 - n_g) (lambda - 1.31) / 1.31, which follows from the definition of the group index and
    reproduces the solved 1290 and 1330 nm indices of the 400 nm strip to 1e-4. `wavelength` may be an
    array; the group index is returned at 1310 nm.
    """
    if not (_STRIP_WIDTHS[0] - 1e-9 <= width <= _STRIP_WIDTHS[-1] + 1e-9):
        raise ValueError(f"strip_index: width {width} um is outside the solved range "
                         f"{_STRIP_WIDTHS[0]} to {_STRIP_WIDTHS[-1]} um")
    n0 = float(np.interp(width, _STRIP_WIDTHS, _STRIP_NEFF))
    g0 = float(np.interp(width, _STRIP_WIDTHS, _STRIP_NGROUP))
    lam = np.asarray(wavelength, float)
    return n0 + (n0 - g0) * (lam - cal.WAVELENGTH) / cal.WAVELENGTH, g0


def _devrec_from_shapes(c, x_range=None):
    """DEVREC from the drawn shapes, with x held to the port plane (see add_drop_ring for why)."""
    boxes = np.array([s.bounds() for shapes in c.get_structures().values() for s in shapes])
    (x0, y0), (x1, y1) = boxes[:, 0].min(axis=0), boxes[:, 1].max(axis=0)
    if x_range is not None:
        x0, x1 = x_range
    c.add("DEVREC", pf.Rectangle((float(x0), float(y0)), (float(x1), float(y1))))


def _extras(pairs):
    """{short name: value} for the parameters that differ from their default, for `_named`."""
    return {key: value for key, value, default in pairs if value != default}


# =====================================================================
# PIN attenuator
# =====================================================================
# [PHYS] Free-carrier absorption at 1310 nm, Nedeljkovic, Soref and Mashanovich, IEEE Photon. J. 3,
# 1171 (2011): d_alpha = 3.48e-22 dN^1.229 + 1.02e-19 dP^1.089, cm^-1 with N and P in cm^-3.
_FCA_E = (3.48e-22, 1.229)
_FCA_H = (1.02e-19, 1.089)
_PIN_TAU = 30e-9           # [SYNTH] effective carrier lifetime in the intrinsic region, s
_PIN_CONFINEMENT = 0.6     # [SYNTH] share of the mode that sits in the carrier-filled region
_PIN_FILM = 0.22           # [PHYS] um, carriers fill the silicon film under the core and slab
_PIN_BASE_LOSS = 0.05      # [SYNTH] dB, implants' tails and the slab tapers at zero current
_PIN_N_MAX = 1.0e20        # [PHYS] cm^-3, the density up to which the Nedeljkovic fit is valid
_Q_E = 1.602176634e-19


def pin_attenuation_db(*, length, intrinsic_width, current):
    """Attenuation in dB that injecting `current` (mA) adds to a `pin_attenuator`. [SYNTH]

    Equal electron and hole density fills the intrinsic region, N = I tau / (q V), with V = length *
    intrinsic_width * 0.22 um and tau = 30 ns, and the loss is the published free-carrier absorption
    at 1310 nm weighted by a confinement of 0.6. The result is a length-scaled loss for one reason:
    the same current in a longer device is a lower density, which is a lower absorption per um, and
    the two nearly cancel in the product, so the loss mostly follows the current, as in a real VOA,
    while a short device pays for it in current density. Roughly 1.2 dB per mA at 500 um by 1.6 um.
    The fit holds up to 1e20 cm^-3, so more current than that raises `ValueError` (it names the limit in mA).
    """
    if current <= 0:
        return 0.0
    volume = length * intrinsic_width * _PIN_FILM * 1e-12          # um^3 -> cm^3
    n = current * 1e-3 * _PIN_TAU / (_Q_E * volume)                # cm^-3
    if n > _PIN_N_MAX:
        raise ValueError(
            f"pin_attenuation_db: {current} mA in {length} um x {intrinsic_width} um injects {n:.1e} cm^-3, above the "
            f"{_PIN_N_MAX:.0e} cm^-3 where the free-carrier-absorption fit holds; use a longer device or less current "
            f"(at most {_PIN_N_MAX * _Q_E * volume / (_PIN_TAU * 1e-3):.1f} mA here)")
    alpha = _PIN_CONFINEMENT * (_FCA_E[0] * n ** _FCA_E[1] + _FCA_H[0] * n ** _FCA_H[1])   # 1/cm, power
    return float(10 / np.log(10) * alpha * length * 1e-4)


@pf.parametric_component(name_prefix="pin")
def pin_attenuator(*, length: pft.PositiveDimension = 200.0,
                   intrinsic_width: pft.PositiveDimension = 1.6,
                   slab_taper: pft.PositiveDimension = 4.0,
                   current: float = 0.0,
                   technology=None):
    """Carrier-injection variable optical attenuator: a rib waveguide with a lateral p-i-n junction.

        P1 <-------------------------------------------------< P0   (rib, 500 nm core on a 90 nm slab)
                .---------------------------------------.
          ANODE |  p+ contact | p | intrinsic | n | n+ contact |  CATHODE     slab taper at both ends
                '---------------------------------------'

    The rib core runs down the middle of an intrinsic gap `intrinsic_width` wide (1.6 um by default);
    p implants fill the slab on the -y side and n implants on the +y side, each with the same graded
    nesting as `pn_phase_shifter` (well, implant, contact implant; IMP.E.1-4), a row of VIA1 on the
    contact implant, and an M1 strip. The ANODE terminal is the p side (y < 0) and CATHODE the n side
    (y > 0), at the input end, as on `pn_phase_shifter`. Forward bias injects carriers into the
    intrinsic region and they absorb light; that is the whole device.

    The slab is 3.0 um wide at both ports, exactly the `SI_RIB_O` spec, and widens over `slab_taper`
    to the width the contacts need (intrinsic_width + 2.3 um), so a rib port plane carries no slab
    steps. P0 is at x = 0 facing +x, P1 at x = `length` facing -x.

    `current` (mA, a cell parameter for the same reason `voltage` is on the modulator: the S matrix is
    cached per component) sets the model: AnalyticWaveguideModel with the rib's index and loss, plus
    `pin_attenuation_db(...)` of free-carrier absorption and 0.05 dB of base loss. It is a trend with
    the right shape and order of magnitude (about 1.2 dB per mA at 500 um by 1.6 um), not a
    prediction, and a current that injects more than 1e20 cm^-3 raises (the absorption fit ends there).
    The attenuation also sits on the cell as `properties["attenuation_db"]`.
    """
    if length < 2 * slab_taper + 6.0:
        raise ValueError(f"pin_attenuator: length {length} um leaves under 6 um between the two slab tapers "
                         f"of {slab_taper} um for contacts")
    if intrinsic_width < 1.0:
        raise ValueError(f"pin_attenuator: intrinsic_width {intrinsic_width} um is narrower than the 1 um "
                         "that keeps the p and n implants clear of the rib's mode")
    tech = _tech(technology)
    spec = _spec(tech, RIB)
    core_w, _ = _profile(spec, tech, "SI")
    slab_w, _ = _profile(spec, tech, "SLAB")
    c = pf.Component(technology=tech)

    ih = intrinsic_width / 2
    sh = ih + 1.15                         # half slab: a 0.95 um contact band plus 0.1 um of nesting each side
    t = slab_taper
    x0, x1 = t, length - t                 # the doped part
    # slab: port width -> contact width -> port width
    c.add("SLAB", pf.Polygon([(0, -slab_w / 2), (t, -sh), (length - t, -sh), (length, -slab_w / 2),
                              (length, slab_w / 2), (length - t, sh), (t, sh), (0, slab_w / 2)]))
    c.add("SI", pf.Rectangle((0, -core_w / 2), (length, core_w / 2)))
    # wells stay apart: the gap between them is the intrinsic region (IMP.O.1 only forbids overlap)
    c.add("PW", pf.Rectangle((x0, -sh), (x1, -ih)), "NW", pf.Rectangle((x0, ih), (x1, sh)))
    c.add("PP", pf.Rectangle((x0, -sh + 0.05), (x1, -ih - 0.05)), "NP", pf.Rectangle((x0, ih + 0.05), (x1, sh - 0.05)))
    c.add("PPP", pf.Rectangle((x0, -sh + 0.10), (x1, -ih - 0.10)), "NPP", pf.Rectangle((x0, ih + 0.10), (x1, sh - 0.10)))
    y_p, y_n = -(sh - 0.575), sh - 0.575   # the middle of each contact band
    n_via = max(1, int((x1 - x0 - 1.0) // 1.0))
    for i in range(n_via):
        x = x0 + 0.5 + i * 1.0
        for y in (y_p, y_n):
            c.add("VIA1", pf.Rectangle(center=(x + 0.15, y), size=(0.3, 0.3)))
    for y in (y_p, y_n):
        c.add("M1", pf.Rectangle((x0 + 0.2, y - 0.35), (x1 - 0.2, y + 0.35)))

    c.add_port(pf.Port((0, 0), 0, spec), "P0")
    c.add_port(pf.Port((length, 0), 180, spec), "P1")
    c.add_terminal({
        "ANODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((x0 + 0.2, y_p - 0.35), (x0 + 4.2, y_p + 0.35))),
        "CATHODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((x0 + 0.2, y_n - 0.35), (x0 + 4.2, y_n + 0.35))),
    })

    attenuation = pin_attenuation_db(length=length, intrinsic_width=intrinsic_width, current=current)
    c.add_model(_wg_model(length, port_spec=RIB, extra_loss=_PIN_BASE_LOSS + attenuation), "analytic")
    c.activate_model("analytic")
    c.properties["attenuation_db"] = attenuation
    _named(c, "pin", l=length, **_extras((("iw", intrinsic_width, 1.6), ("st", slab_taper, 4.0), ("i", current, 0.0))))
    add_recognition(c, icon="eo_ps", devrec=False)
    _devrec_from_shapes(c, (0.0, length))
    return c


# =====================================================================
# Grating antenna
# =====================================================================
_ANT_ALPHA_MAX = 0.02      # [SYNTH] power radiated per um by a duty-0.5 grating of the 70 nm shallow etch, 1/um
_ANT_TIP = 0.16            # [DECK] the tip the PDK's termination cell uses, um (SI.W.1 is 0.15)


def grating_antenna_period(emission_angle=10.0, *, width=1.0, duty=0.75, wavelength=cal.WAVELENGTH):
    """Period (um) that makes a `grating_antenna` of this width and duty radiate at `emission_angle` degrees.

    The grating radiates where the period matches the guided wave to a plane wave in the cladding,
    sin(theta) = n_avg - lambda / period, with n_avg the strip's index averaged over the unetched and
    the etched part of the period (the etched part is 0.295 lower, [SOLVED]).
    """
    n_si, _ = strip_index(width, wavelength)
    n_avg = float(n_si) - (1.0 - duty) * _ETCH_DN
    return round(float(wavelength / (n_avg - np.sin(np.radians(emission_angle)))), 3)   # to the 1 nm grid


@pf.parametric_component(name_prefix="ant")
def grating_antenna(*, period: pft.PositiveDimension = 0.49,
                    duty: float = 0.75,
                    length: pft.PositiveDimension = 100.0,
                    width: pft.PositiveDimension = 1.0,
                    taper_length: pft.PositiveDimension = 10.0,
                    tip_length: pft.PositiveDimension = 10.0,
                    port_spec: str = STRIP, technology=None):
    """Weak grating emitter for a phased array: shallow-etch teeth on a short, wide guide. One port.

        P0 >--taper--[ | | | | | | | | | | | | ]--tip->      (| = 70 nm shallow-etch trench across the guide)

    The guide widens from the strip over `taper_length` to `width`, carries `int(length / period)`
    shallow-etch trenches (the GRAT layer, 70 nm of the 220 nm), and ends in a `tip_length` taper to a
    0.16 um tip, so the end of the guide reflects little. `duty` is the unetched fraction of a period,
    as `fill_factor` is on the grating coupler: the trench is period * (1 - duty) wide. The default
    period 0.49 um with duty 0.75 radiates about 10 degrees off vertical at 1310 nm
    (`grating_antenna_period` solves for any angle). The guide is 1 um wide so that antennas fit a 2 um
    pitch with 1 um between them; a narrower guide radiates a wider beam sideways.

    The trench and the tooth between two trenches must each be 100 nm or more (GRAT.W.1, GRAT.S.1),
    so 0.1 / period <= duty <= 1 - 0.1 / period, and a violation raises here instead of in DRC.

    `width` is limited to 0.4 to 1.0 um, the range the TE0 index (`strip_index`) was solved for. A 1.0 um strip
    also guides TE1 at 1310 nm (cutoff near 0.6 um); the taper from the 0.4 um port keeps the light in TE0, and
    the model is the TE0 one.

    Model: a one-port `TerminationModel` with a reflection of -25 dB, because the cell's job is to
    take power out of the circuit. How much it radiates is stored on the cell, from a simple
    coupled-mode picture: the radiated power coefficient is 0.02 sin(pi duty) per um [SYNTH], so
    `properties["radiated_fraction"]` = 1 - exp(-coefficient * length), `properties["emission_angle_deg"]`
    from the phase matching above, and `properties["divergence_deg"]`, the far-field width along the
    guide, 0.88 lambda / (L_eff cos theta). They are trends, not a prediction.
    """
    trench, tooth = period * (1 - duty), period * duty
    if trench < _RULE_GRAT + _GDS_MARGIN or tooth < _RULE_GRAT + _GDS_MARGIN:
        raise ValueError(
            f"grating_antenna: period {period} um with duty {duty} gives a {trench * 1e3:.0f} nm trench and a "
            f"{tooth * 1e3:.0f} nm tooth; GRAT.W.1 and GRAT.S.1 need 100 nm for both")
    if not 0.4 <= width <= 1.0:
        raise ValueError("grating_antenna: width must be between 0.4 and 1.0 um, the range the TE0 index was solved for")
    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w_port, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    x_g0 = taper_length
    x_g1 = taper_length + length
    x_end = x_g1 + tip_length
    body = pf.Path((0, 0), w_port)
    if abs(width - w_port) > 1e-9:
        body.segment((taper_length, 0), width)
    else:
        body.segment((taper_length, 0))
    body.segment((x_g1, 0)).segment((x_end, 0), _ANT_TIP)
    c.add("SI", body)

    n_teeth = int(length // period)
    for k in range(n_teeth):
        x = x_g0 + k * period + tooth
        c.add("GRAT", pf.Rectangle((x, -width / 2), (x + trench, width / 2)))

    c.add_port(pf.Port((0, 0), 0, spec), "P0")
    c.add_model(pf.TerminationModel(r=10 ** (-25 / 20)), "analytic")
    c.activate_model("analytic")

    n_si, _ = strip_index(width)
    n_avg = float(n_si) - (1 - duty) * _ETCH_DN
    s = n_avg - cal.WAVELENGTH / period
    theta = float(np.degrees(np.arcsin(s))) if abs(s) < 1 else float("nan")
    alpha = _ANT_ALPHA_MAX * np.sin(np.pi * duty)
    l_eff = (1 - np.exp(-alpha * n_teeth * period)) / alpha
    c.properties["emission_angle_deg"] = theta
    c.properties["radiated_fraction"] = float(1 - np.exp(-alpha * n_teeth * period))
    c.properties["divergence_deg"] = float(np.degrees(0.88 * cal.WAVELENGTH / (l_eff * np.cos(np.radians(theta))))) \
        if np.isfinite(theta) else float("nan")
    _named(c, "ant", p=period, d=duty, l=length,
           **_extras((("w", width, 1.0), ("tl", taper_length, 10.0), ("tp", tip_length, 10.0), ("ps", port_spec, STRIP))))
    add_recognition(c, icon="grating_coupler", devrec=False)
    _devrec_from_shapes(c, (0.0, x_end))
    return c


# =====================================================================
# Travelling-wave GSG electrode
# =====================================================================
_TWE_PAD_MARGIN = 2.5      # [DECK] M2 beyond the pad opening, like bond_pad (PAD.E.1 asks 2)
_TWE_VIA = 0.3             # [DECK] VIA3 side (V3.W.1 is 0.25)
# [PHYS] Quasi-static effective permittivity of the line: (eps_Si + 1) / 2, a silicon half space under air. The line
# is 22 um across and only 4.9 um of oxide lies between it and the substrate, so the field is mostly in the silicon;
# the oxide under it lowers this by roughly a tenth, so treat the index as an upper bound (2.5, not the 1.6 of an
# oxide half space) and the impedance (about 50 ohm at 10 / 6 um) as a lower one.
_TWE_EPS_EFF = (11.7 + 1.0) / 2


def _cpw_impedance(signal_width, gap):
    """Quasi-static impedance (ohm) and effective index of an unloaded coplanar line, Hilberg's closed form."""
    k = signal_width / (signal_width + 2 * gap)
    kp = np.sqrt(1 - k * k)
    if k <= np.sqrt(0.5):
        k_over_kp = np.pi / np.log(2 * (1 + np.sqrt(kp)) / (1 - np.sqrt(kp)))
    else:
        k_over_kp = np.log(2 * (1 + np.sqrt(k)) / (1 - np.sqrt(k))) / np.pi
    n_eff = float(np.sqrt(_TWE_EPS_EFF))
    return float(30 * np.pi / n_eff / k_over_kp), n_eff


@pf.parametric_component(name_prefix="twe")
def tw_electrode(*, length: pft.PositiveDimension = 500.0,
                 signal_width: pft.PositiveDimension = 10.0,
                 gap: pft.PositiveDimension = 6.0,
                 ground_width: pft.PositiveDimension = 20.0,
                 arms: int = 2,
                 contact_offset: pft.PositiveDimension = 1.175,
                 tap_pitch: pft.PositiveDimension = 10.0,
                 pad_size: pft.PositiveDimension = 70.0,
                 pad_pitch: pft.PositiveDimension = 100.0,
                 taper_length: pft.PositiveDimension = 100.0,
                 bridge: bool = True,
                 technology=None):
    """Travelling-wave GSG electrode on M2 for a push-pull pair of `pn_phase_shifter` arms, pad sets at both ends.

        G pad                 ground rail (M2)  --------------------------------   G pad
        S pad --taper--  signal rail (M2)       == arm A / arm B in the two gaps ==  --taper-- S pad
        G pad                 ground rail (M2)  --------------------------------   G pad
        x = -taper_length                   x = 0 ........... x = length        x = length + taper_length

    The coplanar line runs along +x from x = 0 to x = `length`: a signal rail `signal_width` wide on
    y = 0 and a ground rail `ground_width` wide at `gap` either side. Each of the two gaps carries one
    modulator arm, centred in it: arm A (y = -(signal_width + gap) / 2) is a `pn_phase_shifter` at
    origin (0, y_A) and arm B another one at (0, +(signal_width + gap) / 2), NOT mirrored. That is the
    push-pull of a dual-arm modulator: arm A has its n side (CATHODE) towards the signal rail and arm B its
    p side (ANODE), so a swing of the signal reverse-biases one arm more and the other less, and the two
    phases move in opposite directions. Mirror arm B (`x_reflection=True`) and both n sides face the rail:
    both arms are reverse-biased together and the phases move the same way, the common-mode drive of two
    phase modulators (`tw_electrode_with_arms(push_pull=False)`). The electrode does not care which way a
    junction faces: the contact strips sit at the same offsets either way. `arms=1` keeps only arm A's taps
    (the single-drive case; the upper ground rail is then only a reference).

    The arms' M1 contact strips sit `contact_offset` (1.175 um, the pn_phase_shifter's) either side of
    each arm's axis. Every `tap_pitch` along the line an M1 jumper 0.7 um wide runs from each strip
    to the rail on its side, ending under it on a VIA3 (the inner strip to the signal rail, the outer
    strip to the ground rail). The rails stay 1.25 um clear of the arm's slab and 3 um clear of its
    core, and `gap` must be at least 3.5 um to keep it so.

    At each end the three rails fan out over `taper_length` to a G-S-G pad set `pad_pitch` apart
    (100 um, a standard probe pitch): M2 squares `pad_size` + 5 um with the PAD opening inside. Four
    ground pads and two signal pads, terminals SIN, GAIN, GBIN (left) and SOUT, GAOUT, GBOUT (right),
    A the -y side and B the +y side. An M1 bridge under the signal taper joins the two ground pads at
    each end (a VIA3 at each side, no via under the signal), so all six ground and signal pads form
    exactly two nets, S and G, as a probe makes them. `bridge=False` leaves the grounds joined only by
    whatever the arms' M1 does.

    There is no DEVREC: like `via_stack` the electrode is wiring, and its footprint lies over the
    arms' devices, which DEV.1 would read as overlapping. The cell has no optical port and no solver
    model; the quasi-static impedance of the unloaded line and its microwave index are stored as
    `properties["z0_ohm"]` and `properties["n_microwave"]` for reference, on the assumption that the field
    sits in the silicon under the line (eps_eff = 6.35, n = 2.52, about 50 ohm at 10 / 6 um): the 4.9 um of
    oxide beneath it lowers both a little, and the arms' capacitance loads the line further, so the loaded
    impedance is lower and the line slower. Do not velocity-match against these numbers; solve the line.
    """
    if arms not in (1, 2):
        raise ValueError("tw_electrode: arms must be 1 or 2")
    if gap < 3.5:
        raise ValueError(f"tw_electrode: gap {gap} um puts the rails closer than 1.25 um to the arm's slab; use 3.5 um or more")
    if length < 20.0:
        raise ValueError("tw_electrode: length under 20 um leaves no room for taps")
    if pad_pitch < pad_size + 2 * _TWE_PAD_MARGIN + 1.0 + 1e-9:
        raise ValueError("tw_electrode: pad_pitch leaves under 1 um between neighbouring pads (M2.S.1)")
    tech = _tech(technology)
    c = pf.Component(technology=tech)
    ws, wg, g, T = signal_width, ground_width, gap, taper_length
    pm = pad_size + 2 * _TWE_PAD_MARGIN
    y_s = ws / 2                       # signal rail edge
    y_g0, y_g1 = ws / 2 + g, ws / 2 + g + wg
    L = length

    # the coplanar line
    c.add("M2", pf.Rectangle((0, -y_s), (L, y_s)),
          pf.Rectangle((0, y_g0), (L, y_g1)), pf.Rectangle((0, -y_g1), (L, -y_g0)))

    # taps: arm A below the signal rail, arm B above it, mirrored so the inner strip is the n side
    ay = (ws + g) / 2
    half = 0.35
    taps = [(-ay, -1)] + ([(+ay, +1)] if arms == 2 else [])
    n_taps = int((L - 10.0) // tap_pitch) + 1
    x_taps = [5.0 + k * tap_pitch for k in range(n_taps)]
    for y_arm, sgn in taps:
        # sgn = -1: the signal rail is at +y of the arm; sgn = +1: at -y of it
        y_in = y_arm - sgn * contact_offset            # strip towards the signal rail
        y_out = y_arm + sgn * contact_offset           # strip towards the ground rail
        y_via_s = sgn * (y_s - 0.5)                    # 0.5 um inside the signal rail's near edge
        y_via_g = sgn * (y_g0 + 0.5)                   # 0.5 um inside the ground rail's near edge
        for x in x_taps:
            for y_from, y_via in ((y_in, y_via_s), (y_out, y_via_g)):
                lo, hi = sorted((y_from, y_via))
                c.add("M1", pf.Rectangle((x - half, lo - half), (x + half, hi + half)))
                c.add("VIA3", pf.Rectangle(center=(x, y_via), size=(_TWE_VIA, _TWE_VIA)))

    # pad sets and the fans that reach them
    pad_ys = {"S": 0.0, "GA": -pad_pitch, "GB": +pad_pitch}
    rail = {"S": (-y_s, y_s), "GA": (-y_g1, -y_g0), "GB": (y_g0, y_g1)}
    terminals = {}
    for end, (x_rail, sign) in {"IN": (0.0, -1), "OUT": (L, +1)}.items():
        x_pad_edge = x_rail + sign * T                 # the pad's edge that meets the fan
        x_pad_c = x_pad_edge + sign * pm / 2
        for name, y_pad in pad_ys.items():
            y_lo, y_hi = rail[name]
            c.add("M2", pf.Polygon([(x_rail, y_lo), (x_rail, y_hi), (x_pad_edge, y_pad + pm / 2), (x_pad_edge, y_pad - pm / 2)]))
            c.add("M2", pf.Rectangle(center=(x_pad_c, y_pad), size=(pm, pm)))
            c.add("PAD", pf.Rectangle(center=(x_pad_c, y_pad), size=(pad_size, pad_size)))
            terminals[f"{name}{end}"] = pf.Terminal(
                _layer(tech, "M2"), pf.Rectangle(center=(x_pad_c, y_pad), size=(pad_size, pad_size)))
        if bridge:
            # M1 under the signal fan, a via on each ground fan: where the fan is half way to the pad
            x_b = x_rail + sign * T / 2
            y_lo_b = 0.5 * (y_g0 + (pad_pitch - pm / 2))
            y_hi_b = 0.5 * (y_g1 + (pad_pitch + pm / 2))
            y_v = 0.5 * (y_lo_b + y_hi_b)
            c.add("M1", pf.Rectangle((x_b - 1.5, -y_v - 1.5), (x_b + 1.5, y_v + 1.5)))
            for yy in (-y_v, y_v):
                for dx in (-0.5, 0.5):
                    for dy in (-0.5, 0.5):
                        c.add("VIA3", pf.Rectangle(center=(x_b + dx, yy + dy), size=(_TWE_VIA, _TWE_VIA)))
    c.add_terminal({name: terminals[name] for name in ("SIN", "GAIN", "GBIN", "SOUT", "GAOUT", "GBOUT")})

    z0, n_mw = _cpw_impedance(ws, g)
    c.properties["z0_ohm"] = z0
    c.properties["n_microwave"] = n_mw
    _named(c, "twe", l=length, ws=signal_width, g=gap, wg=ground_width,
           **_extras((("a", arms, 2), ("co", contact_offset, 1.175), ("tp", tap_pitch, 10.0), ("pz", pad_size, 70.0),
                      ("pp", pad_pitch, 100.0), ("tl", taper_length, 100.0), ("br", bridge, True))))
    return add_recognition(c, icon="tl", devrec=False)


def tw_electrode_with_arms(length=500.0, *, signal_width=10.0, gap=6.0, ground_width=20.0, arms=2, push_pull=True,
                           name=None, **kwargs):
    """A `tw_electrode` placed over its arms: a `pf.Component` holding the electrode and `pn_phase_shifter` arms.

    Arm A is at (0, -(signal_width + gap) / 2) and arm B at (0, +(signal_width + gap) / 2), both `length`
    long, so the arms' optical ports P0 are at x = 0 and P1 at x = `length`, `signal_width + gap` apart.
    With `push_pull=True` (the default) arm B is as it comes, so the signal rail meets the CATHODE of arm A and
    the ANODE of arm B (one arm's reverse bias goes up as the other's goes down); with `push_pull=False` arm B is
    mirrored and the rail meets both CATHODEs (both arms move the same way). The result exposes the electrode's
    six terminals, the arms' terminals as A_ANODE, A_CATHODE, B_ANODE, B_CATHODE, and the arms' ports A0, A1, B0, B1.
    """
    from .cells import pn_phase_shifter

    tech = _tech(kwargs.get("technology"))
    electrode = tw_electrode(length=length, signal_width=signal_width, gap=gap, ground_width=ground_width, arms=arms, **kwargs)
    suffix = "" if (push_pull or arms == 1) else "_cm"
    top = pf.Component(name or f"{electrode.name}_with_arms{suffix}", technology=tech)
    e_ref = top.add_reference(pf.Reference(electrode))
    for term in electrode.terminals:
        top.add_terminal(e_ref[term], term)
    ay = (signal_width + gap) / 2
    arms_ref = {"A": top.add_reference(pf.Reference(pn_phase_shifter(length=length), origin=(0, -ay)))}
    if arms == 2:
        arms_ref["B"] = top.add_reference(pf.Reference(pn_phase_shifter(length=length), origin=(0, ay),
                                                       x_reflection=not push_pull))
    for tag, ref in arms_ref.items():
        top.add_port(ref["P0"], f"{tag}0")
        top.add_port(ref["P1"], f"{tag}1")
        for term in ("ANODE", "CATHODE"):
            top.add_terminal(ref[term], f"{tag}_{term}")
    return top


# =====================================================================
# Contra-directional coupler
# =====================================================================
# Coupled-mode picture of a contra-DC. Guide 1 carries a forward mode (index n1) and guide 2 a
# backward mode (n2); a corrugation of period L couples them when beta1 + beta2 = 2 pi m / L, so
#     L = m lambda / (n1 + n2)
# at the wavelength it drops, and since n1 and n2 differ (the guides have different widths) the
# co-directional coupling between the guides is phase mismatched and does not happen. Light entering
# P0 leaves through P1 except in a stop band of width ~ 2 kappa lambda^2 / (pi (ng1 + ng2)), where it
# leaves at P2 instead, back along the other guide.
_CDC_KAPPA_REF = 0.025     # [SYNTH] 1/um, for a 50 nm tooth, gap 160 nm, order 2, duty 0.32 (the SiEPIC CWDM4 demo's
#                            stop band is a few nm wide at 470 um, which is this order of magnitude)
_CDC_REF = dict(tooth=0.05, gap=0.16, order=2, duty=0.32)
_CDC_BACK = 10 ** (-40 / 20)       # [SYNTH] residual reflection, linear
_CDC_RAD_DB = {1: 0.0, 2: 1.0e-4, 3: 2.0e-4}   # [SYNTH] extra radiation loss of a higher-order grating, dB/um per 50 nm tooth
_CDC_BAND = (1.25, 1.37)           # um, the wavelength span of the stored S matrix
_CDC_STEP = 0.02e-3                # um, 20 pm: a tenth of the side-lobe ripple of a 470 um coupler
_CDC_BLOCK = 40                    # periods per polygon, so no polygon has thousands of vertices


def contra_dc_period(wavelength=1.311, *, order=2, width1=0.40, width2=0.52):
    """Grating period (um, to 0.1 nm) that makes a `contra_dc` drop `wavelength`: m lambda / (n1 + n2).

    `n1` and `n2` are the strip indices of the two widths at `wavelength` (`strip_index`), and `order`
    is the grating order m. The corrugation changes the average width a little, which this leaves out
    (the model of `contra_dc` includes the cell's real ports and tapers, not this shift).
    """
    n1, _ = strip_index(width1, wavelength)
    n2, _ = strip_index(width2, wavelength)
    return round(float(order * wavelength / (n1 + n2)), 4)


def _cdc_kappa(order, duty, gap, tooth):
    """Coupling coefficient (1/um) of a tooth of height `tooth` um: the m-th Fourier part of the tooth,
    the decay of the field across the gap, and a linear dependence on the tooth height. [SYNTH]"""
    c_m = abs(np.sin(np.pi * order * duty)) / (np.pi * order)
    c_ref = abs(np.sin(np.pi * _CDC_REF["order"] * _CDC_REF["duty"])) / (np.pi * _CDC_REF["order"])
    gamma = _field_decay(strip_index(0.46)[0])
    return float(_CDC_KAPPA_REF * (tooth / _CDC_REF["tooth"]) * (c_m / c_ref) * np.exp(-gamma * (gap - _CDC_REF["gap"])))


def _sinhc(sigma, ell):
    """sinh(sigma ell) / sigma, with the sigma -> 0 limit taken exactly."""
    small = np.abs(sigma) < 1e-9
    safe = np.where(small, 1.0, sigma)
    return np.where(small, ell, np.sinh(safe * ell) / safe)


def _cdc_s_array(*, period, order, width1, width2, kappas, end1, end2, lam):
    """4-port S matrix (F, 4, 4), ports P0 P1 (guide 1, left, right) and P2 P3 (guide 2), coupled modes. [SYNTH]

    `kappas` is the coupling of every period; equal values are one uniform grating (a closed form),
    anything else is cut into slices, each solved exactly, and multiplied. `end1` and `end2` are the
    path lengths (um) of the end sections at each end of guide 1 and guide 2.
    """
    kappas = np.asarray(kappas, float)
    n_per = len(kappas)
    lam = np.asarray(lam, float)
    n1, _ = strip_index(width1, lam)
    n2, _ = strip_index(width2, lam)
    alpha = _GUIDE[STRIP][2] * np.log(10) / 10                     # 1/um, power
    rad = _CDC_RAD_DB.get(order, 2e-4 * order / 2) * np.log(10) / 10
    k0 = 2 * np.pi / lam
    b1 = k0 * n1 - 0.5j * (alpha + rad)
    b2 = k0 * n2 - 0.5j * (alpha + rad)
    delta = 0.5 * (b1 + b2) - order * np.pi / period
    ell_total = n_per * period

    m11 = np.ones_like(delta, dtype=complex)
    m12 = np.zeros_like(m11)
    m21 = np.zeros_like(m11)
    m22 = np.ones_like(m11)
    n_slices = 1 if np.ptp(kappas) < 1e-12 else min(n_per, 200)
    for part in np.array_split(np.arange(n_per), n_slices):
        ell, kap = len(part) * period, float(kappas[part].mean())
        sigma = np.lib.scimath.sqrt(kap**2 - delta**2 + 0j)
        c, s = np.cosh(sigma * ell), _sinhc(sigma, ell)
        a11, a12 = c - 1j * delta * s, -1j * kap * s
        a21, a22 = 1j * kap * s, c + 1j * delta * s
        m11, m12, m21, m22 = a11 * m11 + a12 * m21, a11 * m12 + a12 * m22, a21 * m11 + a22 * m21, a21 * m12 + a22 * m22

    t = 1.0 / m22                           # a(L) / a(0) with nothing entering guide 2 at z = L
    r_left = -m21 / m22                     # b(0) / a(0): the light that leaves guide 2 at the input end
    r_right = m12 / m22                     # a(L) / b(L)
    through1 = t * np.exp(1j * (delta - b1) * ell_total)
    through2 = t * np.exp(1j * (delta - b2) * ell_total)
    drop_right = r_right * np.exp(-2j * np.pi * order * ell_total / period)

    s = np.zeros(lam.shape + (4, 4), complex)
    for i in range(4):
        s[..., i, i] = _CDC_BACK
    for i, j, v in ((1, 0, through1), (2, 0, r_left), (1, 3, drop_right), (2, 3, through2)):
        s[..., i, j] = v
        s[..., j, i] = v
    # the end sections: the S-bend and taper at each end, at the mean width's index
    n_e1, _ = strip_index(0.5 * (_STRIP_WIDTHS[0] + width1), lam)
    n_e2, _ = strip_index(0.5 * (_STRIP_WIDTHS[0] + width2), lam)
    e1 = np.exp(-1j * k0 * n_e1 * end1 - 0.5 * alpha * end1)
    e2 = np.exp(-1j * k0 * n_e2 * end2 - 0.5 * alpha * end2)
    e = np.stack([e1, e1, e2, e2], axis=-1)
    return s * e[..., :, None] * e[..., None, :]


def _corrugated_polygons(x0, period, heights, duty, y_lo, y_hi, side, x_stop):
    """Polygons of a straight guide [y_lo, y_hi] whose side `side` (-1 lower, +1 upper) carries one tooth per
    period, `heights[k]` tall, `duty` of the period wide and centred in it, starting at x0 and ending at `x_stop`
    (the grid point the S-bend that follows starts on)."""
    polygons = []
    n = len(heights)
    for b0 in range(0, n, _CDC_BLOCK):
        b1 = min(n, b0 + _CDC_BLOCK)
        xs, xe = x0 + b0 * period, (x_stop if b1 == n else x0 + b1 * period)
        y_e = y_lo if side < 0 else y_hi                  # the corrugated edge
        pts = [(xs, y_hi), (xe, y_hi), (xe, y_lo)] if side < 0 else [(xs, y_lo), (xe, y_lo), (xe, y_hi)]
        for k in range(b1 - 1, b0 - 1, -1):
            t0 = x0 + k * period + (1 - duty) * period / 2
            t1 = t0 + duty * period
            y_t = y_e + side * float(heights[k])
            pts += [(t1, y_e), (t1, y_t), (t0, y_t), (t0, y_e)]
        pts.append((xs, y_e))
        polygons.append(pf.Polygon(pts))
    return polygons


@pf.parametric_component(name_prefix="cdc")
def contra_dc(*, period: pft.PositiveDimension = 0.491,
              num_periods: int = 950,
              width1: pft.PositiveDimension = 0.40,
              width2: pft.PositiveDimension = 0.52,
              gap: pft.PositiveDimension = 0.16,
              corrugation_width: pft.PositiveDimension = 0.05,
              duty: float = 0.32,
              order: int = 2,
              apodization: float = 0.0,
              end_length: pft.PositiveDimension = 20.0,
              port_offset: pft.PositiveDimension = 6.0,
              technology=None):
    """Contra-directional coupler (a Bragg add-drop): two corrugated guides side by side, four ports.

        P0 >--taper--[ guide 1, 400 nm, teeth down ]--taper--> P1        (the through path)
                      [ guide 2, 520 nm, teeth up   ]
        P2 <-s-bend--/                               \\--s-bend-< P3      (guide 2 is bent `port_offset` away)

    Guide 1 (`width1`, 0.40 um, the strip itself) runs straight along y = 0 from P0 at x = 0 to P1.
    Guide 2 (`width2`, 0.52 um, so its index differs by 0.15 and the guides do not couple
    co-directionally) lies under it. Their facing sidewalls carry teeth in phase, `corrugation_width`
    tall, one per `period`, `duty` of the period wide, `gap` apart at the closest (so `gap +
    2 corrugation_width` between the bases). Each end of guide 2 leaves through an S-bend of `end_length`
    that drops it `port_offset` below its axis and tapers it to the 0.4 um port width; guide 1 tapers
    over the same length if `width1` is not 0.4. Ports are the PDK's strip, P0 and P1 on guide 1,
    P2 and P3 on guide 2 (P0, P2 left; P1, P3 right), all facing out. Light in at P0 leaves at P1,
    and in the stop band at P2; light in at P3 leaves at P2, and in the stop band at P1.

    SECOND ORDER, because of the deck. A first-order grating at 1310 nm has a period of 0.246 um, so a
    tooth and a notch are 123 nm, under the 150 nm of SI.W.1 and SI.S.1. At `order` 2 the period doubles
    (`contra_dc_period` gives it) and a tooth of 32 percent of the period is 157 nm, so teeth and notches
    pass; the second harmonic of a pulse train is strongest near 25 percent, so 32 percent still couples
    at 0.9 of the maximum. (`order=1` raises, with the feature size, rather than failing DRC.) A
    second-order grating also radiates a little, which the model carries as 0.1 dB/cm. The period of a
    1271 nm channel is 0.468 um and 0.32 of it is a 149.8 nm tooth, so the short-wavelength channels need a
    duty of 0.33 or more (the CWDM4 build uses 0.34).

    What the model leaves out: each guide also reflects its OWN light, a second-order Bragg grating at
    lambda = n_eff x period (about 1288 nm for guide 1 and 1336 nm for guide 2 at the default period, either
    side of the 1311 nm drop band and inside the stored band; on a 20 nm channel grid they land on the
    neighbouring channels). Only the cross coupling P0-P2 and P1-P3 is modelled, with a flat -40 dB residual
    reflection, so the channel isolation of a chain of these cells is optimistic by those reflections.

    `apodization` (0 to 0.8) shapes the tooth height as (1 - a) + a sin^2(pi z / L), a raised cosine
    that lowers the side lobes; the closest gap stays `gap` at the middle and is larger at the ends.
    `num_periods` x `period` is the grating length (466 um by default).

    Model: `DataModel` of the coupled-mode equations over 1.25 to 1.37 um in 20 pm steps, solved exactly
    slice by slice (one closed form when uniform). Strip indices come from the solved width table
    (`strip_index`), the coupling from `_cdc_kappa` ([SYNTH], 0.025 per um for the default 50 nm tooth),
    loss 2 dB/cm, and the end sections add their phase and loss. Stored on the cell:
    `properties["bragg_wavelength_um"]`, `["kappa_per_um"]` and `["length_um"]`.
    """
    tooth_len, notch_len = duty * period, (1 - duty) * period
    if tooth_len < _RULE_SI + _GDS_MARGIN or notch_len < _RULE_SI + _GDS_MARGIN:
        raise ValueError(
            f"contra_dc: period {period} um with duty {duty} gives a {tooth_len * 1e3:.1f} nm tooth and a "
            f"{notch_len * 1e3:.1f} nm notch; SI.W.1 and SI.S.1 need 153 nm for both (the 150 nm rule plus 3 nm "
            "for the 1 nm grid); use order 2 or more, or a duty nearer the middle")
    if gap < _RULE_SI + _GDS_MARGIN:
        raise ValueError(f"contra_dc: gap {gap * 1e3:.0f} nm is below the 153 nm that SI.S.1 needs (150 nm plus 3 nm for the 1 nm grid)")
    if not 0.0 <= apodization <= 0.8:
        raise ValueError("contra_dc: apodization must be between 0 and 0.8")
    if corrugation_width * (1 - apodization) < 0.02 - 1e-9:
        raise ValueError("contra_dc: the shortest tooth, corrugation_width * (1 - apodization), must be 20 nm or more")
    if abs(width1 - width2) < 0.05:
        raise ValueError("contra_dc: the two guides must differ in width by 50 nm or they couple co-directionally")
    if num_periods < 1:
        raise ValueError("contra_dc: num_periods must be positive")
    tech = _tech(technology)
    spec = _spec(tech, STRIP)
    w_port, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    n_per = int(num_periods)
    grating = n_per * period
    # The grating ends on a grid point (1 nm): the S-bend that follows starts there, and a start that falls
    # between two grid points rounds one way for the S-bend and the other for the last tooth block, which
    # leaves a 1 nm sliver between them, and SI.S.1 reports it.
    x_g0 = end_length
    x_g1 = round(x_g0 + grating, 3)
    x_end = round(x_g1 + end_length, 3)
    u = (np.arange(n_per) + 0.5) / n_per
    heights = corrugation_width * ((1 - apodization) + apodization * np.sin(np.pi * u) ** 2)

    # guide 1 on y = 0, guide 2 below it; both bases leave `gap + 2 corrugation_width` between them
    y1_lo, y1_hi = -width1 / 2, width1 / 2
    y2_hi = y1_lo - (gap + 2 * corrugation_width)
    y2_lo = y2_hi - width2
    y2_c = 0.5 * (y2_lo + y2_hi)
    for poly in _corrugated_polygons(x_g0, period, heights, duty, y1_lo, y1_hi, -1, x_g1):
        c.add("SI", poly)
    for poly in _corrugated_polygons(x_g0, period, heights, duty, y2_lo, y2_hi, +1, x_g1):
        c.add("SI", poly)

    # ends: guide 1 straight (tapered to the port width if it is not 0.4 um), guide 2 bent away
    for x_a, x_b in ((0.0, x_g0), (x_end, x_g1)):
        c.add("SI", pf.Path((x_a, 0.0), w_port).segment((x_b, 0.0), width1))
    left = pf.Path((0.0, y2_c - port_offset), w_port).s_bend((end_length, port_offset), 0.2, direction=(1, 0), width=width2, relative=True)
    right = pf.Path((x_g1, y2_c), width2).s_bend((end_length, -port_offset), 0.2, direction=(1, 0), width=w_port, relative=True)
    c.add("SI", left, right)
    s_len = float(left.length())

    c.add_port(pf.Port((0.0, 0.0), 0, spec), "P0")
    c.add_port(pf.Port((x_end, 0.0), 180, spec), "P1")
    c.add_port(pf.Port((0.0, y2_c - port_offset), 0, spec), "P2")
    c.add_port(pf.Port((x_end, y2_c - port_offset), 180, spec), "P3")

    kappas = np.array([_cdc_kappa(order, duty, gap, h) for h in heights])
    lam = np.arange(_CDC_BAND[0], _CDC_BAND[1] + 0.5 * _CDC_STEP, _CDC_STEP)
    s_array = _cdc_s_array(period=period, order=order, width1=width1, width2=width2, kappas=kappas,
                           end1=end_length, end2=s_len, lam=lam)
    freqs = pf.C_0 / lam
    idx = np.argsort(freqs)
    c.add_model(pf.DataModel(s_array=s_array[idx], frequencies=freqs[idx], ports=["P0", "P1", "P2", "P3"]), "coupled_mode")
    c.activate_model("coupled_mode")

    # the wavelength where m lambda = period (n1 + n2). Both indices are linear in lambda, so this is exact:
    # lambda_B = period (s0 - slope lambda0) / (m - period slope), s0 and slope the sum of the two indices and its slope at 1310 nm
    n1_0, g1 = strip_index(width1)
    n2_0, g2 = strip_index(width2)
    slope = (float(n1_0 + n2_0) - (g1 + g2)) / cal.WAVELENGTH
    lam_b = period * (float(n1_0 + n2_0) - slope * cal.WAVELENGTH) / (order - period * slope)
    c.properties["bragg_wavelength_um"] = float(lam_b)
    c.properties["kappa_per_um"] = float(kappas.max())
    c.properties["length_um"] = float(grating)
    _named(c, "cdc", p=period, n=n_per, w2=width2,
           **_extras((("w1", width1, 0.40), ("g", gap, 0.16), ("dw", corrugation_width, 0.05), ("d", duty, 0.32),
                      ("m", order, 2), ("a", apodization, 0.0), ("el", end_length, 20.0), ("po", port_offset, 6.0))))
    add_recognition(c, icon="wdm", devrec=False)
    _devrec_from_shapes(c, (0.0, x_end))
    return c
