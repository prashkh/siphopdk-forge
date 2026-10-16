"""SiPhoPDK Tier-1 device library.

Every cell is a `@pf.parametric_component`, draws the recognition layers, carries a model
that runs offline (analytic, or a DataModel built from synthetic measurement data), and is
clean against `rules/siphopdk.drc`. No cell requires a simulation to be usable.

Cells that exist in the parametric library are wrapped rather than redrawn, so the PDK
inherits the Euler bends, path bookkeeping and port handling that already work.
"""

import numpy as np
import photonforge as pf
import photonforge.typing as pft

from . import _calibration as cal
from ._recognition import add_recognition
from .technology import siphopdk as _tech_fn

STRIP = "SI_STRIP_O"
RIB = "SI_RIB_O"
WIDE = "SI_WIDE_O"
MM = "SI_MM_O"
SIN = "SIN_STRIP_O"

# Where the ring modulator's heater leads leave the ring. Symmetric about the vertical so
# the top of the device stays free for the cathode escape, and far enough apart that the
# two leads reach two pads without either of them crossing the other's route.
_LEAD_ANGLES = (135.0, 45.0)


def _tech(technology=None):
    """The technology a SiPhoPDK cell is built on.

    `pf.config.default_technology` is pre-populated with a "Basic Technology", so trusting
    it blindly builds PDK cells on the wrong layer table and fails with a bare KeyError on
    the first port-spec lookup. Only honour the default when it really is this PDK.
    """
    if technology is not None:
        return technology
    default = pf.config.default_technology
    if default is not None and default.name == "SiPhoPDK":
        return default
    return _tech_fn()


def _spec(technology, name):
    return _tech(technology).ports[name]


def _profile(spec, technology, layer):
    """A port spec's (width, offset) on a layer of this cell's technology, named.

    `_profile(PortSpec, tech, "SI")` looks the name up in photonforge's DEFAULT technology,
    not the cell's: a cell built with `technology=` while the default was something else
    (the LiveViewer's component library builds PDK cells that way) failed with "Layer name
    not found in technology". The layer is resolved here, through the cell's own technology.
    """
    return spec.path_profile_for(_layer(technology, layer), technology)


def _layer(technology, name):
    """A layer of this cell's technology as (layer, datatype). `pf.Terminal(_layer(tech, "M1"), ...)` and
    `path_profile_for("SI")` resolve a name through the DEFAULT technology instead."""
    return tuple(technology.layers[name].layer)


def _named(c, prefix, **params):
    """Readable palette name: the PDK prefix plus the parameters that matter.

    A cell built by wrapping `pf.parametric.*` otherwise keeps that function's own hashed
    name, so the palette shows `s_bend_straight_coupler__Y4ZJ...` instead of `dc_...`.
    """
    bits = "_".join(f"{k}{v}".replace(".", "p") for k, v in params.items())
    c.name = f"{prefix}_{bits}" if bits else prefix
    return c


# The index and loss of each guiding cross-section, keyed by the port spec that names it.
# Written as a table rather than a `rib=True` flag because the flag only had two states and
# the PDK has three guides: `straight(port_spec="SIN_STRIP_O")` was accepted, drawn correctly
# on the nitride layer, and then modelled with silicon's index - a SiN waveguide reporting
# n_eff 2.79 instead of 1.64, which is a 70 percent error in delay and nothing anywhere
# complains. Every cell that takes a `port_spec` now gets its model from this.
_GUIDE = {
    STRIP: (cal.N_EFF_STRIP, cal.N_GROUP_STRIP, cal.LOSS_STRIP),
    RIB: (cal.N_EFF_RIB, cal.N_GROUP_RIB, cal.LOSS_RIB),
    WIDE: (cal.N_EFF_STRIP, cal.N_GROUP_STRIP, cal.LOSS_STRIP),
    MM: (cal.N_EFF_STRIP, cal.N_GROUP_STRIP, cal.LOSS_STRIP),
    SIN: (cal.N_EFF_SIN, cal.N_GROUP_SIN, cal.LOSS_SIN),
}


def _wg_model(length, *, port_spec=STRIP, extra_loss=0.0, **kwargs):
    try:
        n_eff, n_g, loss = _GUIDE[port_spec]
    except KeyError:
        raise KeyError(
            f"no guide constants for port spec {port_spec!r}; add it to cells._GUIDE rather "
            "than letting it default, or the cell is drawn on one layer and modelled as "
            "another"
        ) from None
    return pf.AnalyticWaveguideModel(
        n_eff=n_eff, n_group=n_g, length=length, propagation_loss=loss,
        extra_loss=extra_loss, reference_frequency=cal.REF_FREQUENCY, **kwargs
    )


# =====================================================================
# Passives
# =====================================================================
@pf.parametric_component(name_prefix="straight")
def straight(*, length: pft.PositiveDimension = 10.0, port_spec: str = STRIP, technology=None):
    """Straight waveguide."""
    tech = _tech(technology)
    c = pf.parametric.straight(port_spec=_spec(tech, port_spec), length=length, technology=tech)
    c.add_model(_wg_model(length, port_spec=port_spec), "analytic")
    c.activate_model("analytic")
    _named(c, "straight", l=length)
    return add_recognition(c, icon="wg")


@pf.parametric_component(name_prefix="bend")
def bend(*, radius: pft.PositiveDimension = 10.0, angle: pft.Angle = 90.0,
         port_spec: str = STRIP, euler_fraction: float = 0.2, technology=None):
    """Euler-assisted circular bend."""
    tech = _tech(technology)
    c = pf.parametric.bend(port_spec=_spec(tech, port_spec), radius=radius, angle=angle,
                           euler_fraction=euler_fraction, technology=tech)
    arc = abs(np.deg2rad(angle)) * radius
    c.add_model(_wg_model(arc, port_spec=port_spec,
                          extra_loss=cal.LOSS_BEND_90 * abs(angle) / 90), "analytic")
    c.activate_model("analytic")
    _named(c, "bend", r=radius, a=angle)
    return add_recognition(c, icon="bend")


@pf.parametric_component(name_prefix="taper")
def taper(*, length: pft.PositiveDimension = 20.0, port_spec1: str = STRIP,
          port_spec2: str = WIDE, technology=None):
    """Linear width taper between two port specs on the same layer.

    The default target is the 1 um wide routing waveguide, which the PDK declares with a
    single mode: an analytic waveguide model needs the same mode count at both ports, and a
    taper feeding a multimode spec has to carry a multimode model instead.
    """
    tech = _tech(technology)
    c = pf.parametric.transition(port_spec1=_spec(tech, port_spec1),
                                 port_spec2=_spec(tech, port_spec2),
                                 length=length, technology=tech)
    c.add_model(_wg_model(length, extra_loss=0.01), "analytic")
    c.activate_model("analytic")
    _named(c, "taper", l=length)
    return add_recognition(c, icon="taper")


@pf.parametric_component(name_prefix="s2r")
def strip_to_rib(*, length: pft.PositiveDimension = 20.0,
                 slab_launch: pft.PositiveDimension = 1.5,
                 slab_start: pft.PositiveDimension = 1.0, technology=None):
    """Strip-to-rib transition.

    The modulator and the detector are rib devices, the routing layer is a strip, so every
    rib port needs one of these. Two details are set by the deck rather than by taste: the
    slab starts a little inside the cell so the strip port's cross-section really is a bare
    strip (a port spec has to list every mask layer present at the port), and it launches at
    `slab_launch` wide so it already encloses the core by more than 500 nm and never trips
    SLAB.E.1.
    """
    tech = _tech(technology)
    strip, rib = _spec(tech, STRIP), _spec(tech, RIB)
    w0, _ = _profile(strip, tech, "SI")
    w1, _ = _profile(rib, tech, "SI")
    slab_w, _ = _profile(rib, tech, "SLAB")
    c = pf.Component(technology=tech)
    c.add("SI", pf.Path((0, 0), w0).segment((length, 0), w1))
    c.add("SLAB", pf.Path((slab_start, 0), slab_launch).segment((length, 0), slab_w))
    c.add_port(pf.Port((0, 0), 0, strip), "P0")
    c.add_port(pf.Port((length, 0), 180, rib), "P1")
    c.add_model(_wg_model(length, extra_loss=0.02), "analytic")
    c.activate_model("analytic")
    _named(c, "s2r", l=length)
    return add_recognition(c, icon="transition")


@pf.parametric_component(name_prefix="crossing")
def crossing(*, arm_length: pft.PositiveDimension = 6.0, added_width: pft.Dimension = 0.8,
             port_spec: str = STRIP, technology=None):
    """Waveguide crossing built from two tapered arms."""
    tech = _tech(technology)
    c = pf.parametric.crossing(port_spec=_spec(tech, port_spec), arm_length=arm_length,
                               added_width=added_width, technology=tech)
    t = 10 ** (-cal.LOSS_CROSSING / 20)
    c.add_model(pf.CrossingModel(t=t, x=cal.CROSSTALK_CROSSING, r=1e-3,
                                 propagation_length=2 * arm_length,
                                 n_eff=cal.N_EFF_STRIP, n_group=cal.N_GROUP_STRIP,
                                 reference_frequency=cal.REF_FREQUENCY), "analytic")
    c.activate_model("analytic")
    _named(c, "crossing", a=arm_length)
    return add_recognition(c, icon="crossing")


@pf.parametric_component(name_prefix="ybranch")
def y_branch(*, body_length: pft.PositiveDimension = 3.0, body_width: pft.PositiveDimension = 1.4,
             arm_pitch: pft.PositiveDimension = 0.7, s_bend_length: pft.PositiveDimension = 12.0,
             arm_offset: pft.PositiveDimension = 2.5, port_spec: str = STRIP, technology=None):
    """1x2 splitter: a tapered body that hands off to two s-bends.

    The arms leave the body already `arm_pitch` apart, so the gap between them never passes
    through the sub-150 nm range that would break SI.S.1. That is a deck-driven design
    choice, which is the point of authoring the deck first.
    """
    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    body = (pf.Path((-2.0, 0), w).segment((0, 0)).segment((body_length, 0), body_width))
    arms = []
    for sign in (1, -1):
        arms.append(
            pf.Path((body_length, sign * arm_pitch / 2), w)
            .s_bend((s_bend_length, sign * (arm_offset - arm_pitch / 2)), 0.2, relative=True)
            .segment((2.0, 0), relative=True)
        )
    c.add("SI", body, *arms)
    x_out = body_length + s_bend_length + 2.0
    c.add_port(pf.Port((-2.0, 0), 0, spec), "P0")
    c.add_port(pf.Port((x_out, arm_offset - arm_pitch / 2 + arm_pitch / 2), 180, spec), "P1")
    c.add_port(pf.Port((x_out, -(arm_offset - arm_pitch / 2) - arm_pitch / 2), 180, spec), "P2")
    t = 10 ** (-0.15 / 20) / np.sqrt(2)  # [SYNTH] 0.15 dB excess loss
    c.add_model(pf.PowerSplitterModel(t=t, r0=1e-3, r1=1e-3), "analytic")
    c.activate_model("analytic")
    _named(c, "ybranch")
    return add_recognition(c, icon="y-splitter")


@pf.parametric_component(name_prefix="mmi2x2")
def mmi_2x2(*, body_length: pft.PositiveDimension = 30.0, body_width: pft.PositiveDimension = 3.5,
            access_pitch: pft.PositiveDimension = 1.8, access_width: pft.PositiveDimension = 1.0,
            taper_length: pft.PositiveDimension = 5.0, port_spec: str = STRIP, technology=None):
    """2x2 multimode interference coupler, nominally 50:50."""
    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    body = pf.Rectangle((0, -body_width / 2), (body_length, body_width / 2))
    access = []
    for sign in (1, -1):
        y = sign * access_pitch / 2
        access.append(pf.Path((-taper_length, y), w).segment((0, y), access_width))
        access.append(pf.Path((body_length + taper_length, y), w).segment((body_length, y), access_width))
    c.add("SI", body, *access)
    for i, (x, d, sign) in enumerate(
        ((-taper_length, 0, 1), (-taper_length, 0, -1),
         (body_length + taper_length, 180, 1), (body_length + taper_length, 180, -1))
    ):
        c.add_port(pf.Port((x, sign * access_pitch / 2), d, spec), f"P{i}")
    # [SYNTH] 50:50 with 0.3 dB excess loss and the usual -90 deg cross phase
    a = 10 ** (-0.3 / 20) / np.sqrt(2)
    c.add_model(pf.DirectionalCouplerModel(t=a, c=a * 1j, i=1e-3, r=1e-3,
                                           propagation_length=body_length + 2 * taper_length,
                                           n_eff=cal.N_EFF_STRIP, n_group=cal.N_GROUP_STRIP,
                                           reference_frequency=cal.REF_FREQUENCY), "analytic")
    c.activate_model("analytic")
    _named(c, "mmi2x2", l=body_length)
    return add_recognition(c, icon="mmi2x2")


@pf.parametric_component(name_prefix="dc")
def directional_coupler(*, coupling_distance: pft.PositiveDimension = 0.6,
                        coupling_length: pft.PositiveDimension = 15.0,
                        s_bend_length: pft.PositiveDimension = 15.0,
                        s_bend_offset: pft.PositiveDimension = 4.0,
                        port_spec: str = STRIP, technology=None):
    """Symmetric directional coupler: both arms fan out through mirrored s-bends.

    `s_bend_coupler` bends both sides; `s_bend_straight_coupler` leaves one arm straight,
    which makes the device asymmetric and the two outputs physically inequivalent.
    """
    tech = _tech(technology)
    c = pf.parametric.s_bend_coupler(
        port_spec=_spec(tech, port_spec), coupling_distance=coupling_distance,
        coupling_length=coupling_length, s_bend_length=s_bend_length,
        s_bend_offset=s_bend_offset, technology=tech)
    # [SYNTH] coupling length for full power transfer at the design gap
    c.add_model(pf.AnalyticDirectionalCouplerModel(
        interaction_length=coupling_length, coupling_length=30.0,
        propagation_length=coupling_length + 2 * s_bend_length,
        # insertion_loss is in dB, but isolation and reflection are linear
        # coefficients: the same constructor mixes both conventions.
        insertion_loss=0.1, isolation=10 ** (-40 / 20), reflection=10 ** (-60 / 20),
        n_eff=cal.N_EFF_STRIP, n_group=cal.N_GROUP_STRIP,
        reference_frequency=cal.REF_FREQUENCY), "analytic")
    c.activate_model("analytic")
    _named(c, "dc", lc=coupling_length, g=coupling_distance)
    return add_recognition(c, icon="dc")


def _gc_data_model(port_names):
    """Synthetic measurement data for the grating coupler.

    Shaped like the SPGC measured in the platform paper: Gaussian in wavelength, 1.3 dB peak
    loss, 20 nm 1 dB bandwidth. This stands in for a measured .s2p until real data exists,
    and it exercises the same DataModel path a measured device would use.
    """
    lam = np.linspace(1.26, 1.36, 201)
    sigma = cal.GC_BW_1DB_NM * 1e-3 / (2 * np.sqrt(2 * np.log(10 ** (1 / 10))))
    loss_db = cal.GC_PEAK_LOSS_DB + 10 * np.log10(np.e) * ((lam - cal.GC_PEAK_WAVELENGTH) ** 2) / (2 * sigma**2)
    t = 10 ** (-loss_db / 20)
    r = 10 ** (-25 / 20)
    freqs = pf.C_0 / lam
    order = np.argsort(freqs)
    n = lam.size
    s_array = np.zeros((n, 2, 2), dtype=complex)
    s_array[:, 0, 0] = r
    s_array[:, 1, 1] = r
    s_array[:, 0, 1] = t[order]
    s_array[:, 1, 0] = t[order]
    return pf.DataModel(s_array=s_array, frequencies=freqs[order], ports=list(port_names))


@pf.parametric_component(name_prefix="gc")
def grating_coupler(*, taper_length: pft.PositiveDimension = 25.0,
                    aperture: pft.PositiveDimension = 14.0,
                    period: pft.PositiveDimension = 0.58, fill_factor: float = 0.5,
                    num_teeth: int = 26, stub_length: pft.PositiveDimension = 2.0,
                    port_spec: str = STRIP, technology=None):
    """Focusing 1D grating coupler, O-band, with a synthetic measured model.

    The teeth are arcs concentric with the waveguide tip rather than straight lines, so the
    fan itself focuses the beam and no separate slab taper is needed. The arc centre sits a
    distance `w / (2 tan(theta))` behind the waveguide end, which is where the sector chord
    equals the waveguide width, so the fan grows out of the waveguide instead of starting
    from a sliver.

    `aperture` is the chord at the first tooth; the fan keeps opening across the grating, so
    the far end is wider. The default 14 um over a 25 um taper gives a 15.6 degree half
    angle, which is the usual proportion for a focusing coupler to a 10.4 um fibre mode.
    """
    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    half_angle = np.degrees(np.arctan((aperture / 2) / taper_length))
    d = w / (2 * np.tan(np.radians(half_angle)))  # arc centre, behind the waveguide tip
    centre = (-d, 0.0)
    sector = [-half_angle, half_angle]
    r_first = taper_length
    r_last = r_first + num_teeth * period
    r_out = r_last + period

    fan = pf.Circle(radius=r_out, center=centre, sector=sector)
    # the stub covers the narrow part of the fan, up to where the chord reaches w
    stub = pf.Rectangle((-d - stub_length, -w / 2), (0.5, w / 2))
    c.add("SI", fan, stub)

    # Each tooth is intersected with the fan: pf.Circle approximates an arc with chords, and
    # the chord counts at two different radii do not nest, so a raw annular sector pokes
    # outside the silicon at the sector corners and trips GRAT.I.1.
    trench = period * (1 - fill_factor)
    for i in range(num_teeth):
        r = r_first + i * period
        tooth = pf.Circle(radius=r + trench, inner_radius=r, center=centre, sector=sector)
        c.add("GRAT", *pf.boolean(tooth, fan, "*"))

    c.add_port(pf.Port((-d - stub_length, 0), 0, spec), "P0")
    x_fibre = -d + 0.5 * (r_first + r_last)
    c.add_port(pf.FiberPort((x_fibre, 0, 8.0), (0.17, 0, -1), (10.4, 10.4), (-4.0, 4.0)), "FIBER")
    c.add_model(_gc_data_model(("P0", "FIBER")), "measured")
    c.activate_model("measured")
    _named(c, "gc", n=num_teeth, p=period)
    return add_recognition(c, icon="grating_coupler")


@pf.parametric_component(name_prefix="term")
def termination(*, length: pft.PositiveDimension = 25.0,
                tip_width: pft.PositiveDimension = 0.16,
                port_spec: str = STRIP, technology=None):
    """Radiative termination: a slow taper down to a tip that is still printable.

    Every unused arm of a coupler needs one, otherwise the extracted netlist carries a
    dangling port and the reflection is unmodelled.
    """
    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)
    c.add("SI", pf.Path((0, 0), w).segment((length, 0), tip_width))
    c.add_port(pf.Port((0, 0), 0, spec), "P0")
    c.add_model(pf.TerminationModel(r=10 ** (-35 / 20)), "analytic")
    c.activate_model("analytic")
    _named(c, "term", l=length)
    return add_recognition(c, icon="termination")


@pf.parametric_component(name_prefix="tap")
def tap(*, coupling_distance: pft.PositiveDimension = 0.6,
        coupling_length: pft.PositiveDimension = 2.0,
        s_bend_length: pft.PositiveDimension = 25.0,
        s_bend_offset: pft.PositiveDimension = 12.0,
        port_spec: str = STRIP, technology=None):
    """Weak monitor tap, nominally 1 % to the drop arm.

    Deliberately asymmetric: `s_bend_straight_coupler` keeps the through line straight, which
    is what a tap wants. The symmetric `s_bend_coupler` is the right choice for a 50:50
    splitter, not for a monitor.

    The arms separate by 12 um rather than the few microns the coupling needs, because the
    unused arm gets a butt-connected terminator and at 4 um that terminator sits inside the
    router's approach corridor to the input port.
    """
    tech = _tech(technology)
    c = pf.parametric.s_bend_straight_coupler(
        port_spec=_spec(tech, port_spec), coupling_distance=coupling_distance,
        coupling_length=coupling_length, s_bend_length=s_bend_length,
        s_bend_offset=s_bend_offset, technology=tech)
    # cross power = sin^2(pi/2 * l_i / l_c); l_i = 2 um against a 30 um beat length gives 1.1 %
    c.add_model(pf.AnalyticDirectionalCouplerModel(
        interaction_length=coupling_length, coupling_length=30.0,
        propagation_length=coupling_length + 2 * s_bend_length,
        insertion_loss=0.05, isolation=10 ** (-40 / 20), reflection=10 ** (-60 / 20),
        n_eff=cal.N_EFF_STRIP, n_group=cal.N_GROUP_STRIP,
        reference_frequency=cal.REF_FREQUENCY), "analytic")
    c.activate_model("analytic")
    _named(c, "tap", lc=coupling_length, g=coupling_distance)
    return add_recognition(c, icon="dc")


@pf.parametric_component(name_prefix="pd")
def photodiode(*, length: pft.PositiveDimension = 25.0,
               ge_width: pft.PositiveDimension = 2.0,
               slab_half_width: pft.PositiveDimension = 4.0,
               input_length: pft.PositiveDimension = 5.0,
               technology=None):
    """Evanescently coupled germanium photodiode, lateral p-i-n on the slab.

    The rib core runs under the germanium; the p and n contacts sit on the slab on either
    side, so no contact lands on the germanium itself. Responsivity, dark current and
    bandwidth follow the platform paper.
    """
    tech = _tech(technology)
    spec = _spec(tech, RIB)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)
    x0, x1 = -input_length, length
    sh = slab_half_width

    c.add("SLAB", pf.Rectangle((x0, -sh), (x1, sh)))
    c.add("SI", pf.Rectangle((x0, -w / 2), (x1 - 1.0, w / 2)))
    c.add("GE", pf.Rectangle((2.0, -ge_width / 2), (x1 - 2.0, ge_width / 2)))
    # wells kept apart (IMP.O.1), each graded implant inside the previous one (IMP.E.1-4)
    c.add("PW", pf.Rectangle((x0, -sh), (x1, -0.6)), "NW", pf.Rectangle((x0, 0.6), (x1, sh)))
    c.add("PP", pf.Rectangle((x0, -sh + 0.1), (x1, -0.7)),
          "NP", pf.Rectangle((x0, 0.7), (x1, sh - 0.1)))
    c.add("PPP", pf.Rectangle((x0, -sh + 0.2), (x1, -sh + 1.1)),
          "NPP", pf.Rectangle((x0, sh - 1.1), (x1, sh - 0.2)))
    y_p, y_n = -sh + 0.65, sh - 0.65
    n_via = max(1, int((x1 - 1.0) // 1.0))
    for i in range(n_via):
        x = 0.5 + i * 1.0
        for y in (y_p, y_n):
            c.add("VIA1", pf.Rectangle(center=(x, y), size=(0.3, 0.3)))
    for y in (y_p, y_n):
        c.add("M1", pf.Rectangle((0.2, y - 0.35), (x1 - 0.2, y + 0.35)))

    c.add_port(pf.Port((x0, 0), 0, spec), "P0")
    c.add_port(pf.Port((0.5 * (x0 + x1), sh + 1.0), -90,
                       pf.virtual_port_spec(classification="electrical")), "OUT")
    c.add_terminal({
        "ANODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((0.2, y_p - 0.35), (4.0, y_p + 0.35))),
        "CATHODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((0.2, y_n - 0.35), (4.0, y_n + 0.35))),
    })
    model = pf.TerminationModel(r=10 ** (-30 / 20))
    model.time_stepper = pf.PhotodiodeTimeStepper(
        responsivity=cal.PD_RESPONSIVITY, dark_current=cal.PD_DARK_CURRENT,
        filter_frequency=cal.PD_BANDWIDTH, gain=1.0, z0=50.0)
    c.add_model(model, "photodiode")
    c.activate_model("photodiode")
    c.add_model(pf.TerminationModel(r=0), "electrical", set_active="electrical")
    _named(c, "pd", l=length)
    return add_recognition(c, icon="photodiode")


@pf.parametric_component(name_prefix="via")
def via_stack(*, size: pft.PositiveDimension = 2.0,
              m1_margin: pft.PositiveDimension = 0.4,
              m2_margin: pft.PositiveDimension = 0.5, technology=None):
    """M1-to-M2 via stack.

    Signals run on M1 from the devices and the ground rail runs on M2 above them, so every
    ground tap needs one of these. Margins come straight from the deck: M1 must enclose
    VIA3 by 100 nm and M2 by 150 nm, and both are given more than that.
    """
    tech = _tech(technology)
    c = pf.Component(technology=tech)
    h = size / 2
    c.add("VIA3", pf.Rectangle((-h, -h), (h, h)))
    c.add("M1", pf.Rectangle((-h - m1_margin, -h - m1_margin), (h + m1_margin, h + m1_margin)))
    c.add("M2", pf.Rectangle((-h - m2_margin, -h - m2_margin), (h + m2_margin, h + m2_margin)))
    c.add_terminal({
        "M1": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((-h - m1_margin, -h - m1_margin),
                                             (h + m1_margin, h + m1_margin))),
        "M2": pf.Terminal(_layer(tech, "M2"), pf.Rectangle((-h - m2_margin, -h - m2_margin),
                                             (h + m2_margin, h + m2_margin))),
    })
    _named(c, "via", s=size)
    return add_recognition(c, icon="connection", devrec=False)


@pf.parametric_component(name_prefix="pad")
def bond_pad(*, pad_size: pft.PositiveDimension = 80.0,
             m2_margin: pft.PositiveDimension = 2.5,
             via_size: pft.PositiveDimension = 2.0, via_pitch: pft.PositiveDimension = 5.0,
             technology=None):
    """Bond pad: a top-metal square with a passivation opening, on an M1 landing.

    The M1 landing and the via array underneath let a signal trace arrive on M1, which is
    what the device terminals use, so a route does not need its own via stack at the pad end.
    Two terminals are exposed, one per metal.
    """
    tech = _tech(technology)
    c = pf.Component(technology=tech)
    h = pad_size / 2
    hl = h - 2.0
    c.add("M2", pf.Rectangle((-h - m2_margin, -h - m2_margin), (h + m2_margin, h + m2_margin)))
    c.add("PAD", pf.Rectangle((-h, -h), (h, h)))
    c.add("M1", pf.Rectangle((-hl, -hl), (hl, hl)))
    n = max(1, int((2 * hl - via_size) // via_pitch))
    x0 = -0.5 * (n - 1) * via_pitch
    for i in range(n):
        for j in range(n):
            c.add("VIA3", pf.Rectangle(center=(x0 + i * via_pitch, x0 + j * via_pitch),
                                       size=(via_size, via_size)))
    c.add_terminal({
        "T": pf.Terminal(_layer(tech, "M2"), pf.Rectangle((-h, -h), (h, h))),
        "M1": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((-hl, -hl), (hl, hl))),
    })
    _named(c, "pad", s=pad_size)
    return add_recognition(c, icon="bondpad")


# =====================================================================
# Actives
# =====================================================================
@pf.parametric_component(name_prefix="tps")
def thermal_phase_shifter(*, length: pft.PositiveDimension = 100.0,
                          heater_width: pft.PositiveDimension = 2.0,
                          heater_offset: pft.Dimension = 0.0,
                          pad_size: pft.PositiveDimension = 4.0,
                          port_spec: str = STRIP, technology=None):
    """Strip waveguide with a TiN heater above it and M1 landing pads."""
    tech = _tech(technology)
    spec = _spec(tech, port_spec)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)

    c.add("SI", pf.Rectangle((0, -w / 2), (length, w / 2)))
    y0 = heater_offset - heater_width / 2
    c.add("HTR", pf.Rectangle((0, y0), (length, y0 + heater_width)))
    # a via + M1 pad at each end of the heater
    for x in (0.0, length - pad_size):
        c.add("HTR", pf.Rectangle((x, y0 - pad_size), (x + pad_size, y0 + heater_width + pad_size)))
        c.add("VIA2", pf.Rectangle(center=(x + pad_size / 2, heater_offset), size=(0.4, 0.4)))
        c.add("M1", pf.Rectangle(center=(x + pad_size / 2, heater_offset),
                                 size=(pad_size - 0.5, pad_size - 0.5)))
    c.add_port(pf.Port((0, 0), 0, spec), "P0")
    c.add_port(pf.Port((length, 0), 180, spec), "P1")
    c.add_terminal({
        "HTRP": pf.Terminal(_layer(tech, "M1"), pf.Rectangle(center=(pad_size / 2, heater_offset),
                                             size=(pad_size - 0.5, pad_size - 0.5))),
        "HTRN": pf.Terminal(_layer(tech, "M1"), pf.Rectangle(center=(length - pad_size / 2, heater_offset),
                                             size=(pad_size - 0.5, pad_size - 0.5))),
    })
    c.add_model(_wg_model(length, dn_dT=cal.DN_DT, temperature=293.0,
                          reference_temperature=293.0), "analytic")
    c.activate_model("analytic")
    _named(c, "tps", l=length)
    return add_recognition(c, icon="to_ps")


def _junction_rows(c, length, *, slab_half=1.75, core_half=0.25, via_pitch=1.0):
    """Draw a lateral p-n junction along a rib: wells, graded implants, contacts, M1.

    Geometry is chosen against the deck: wells never overlap (IMP.O.1), each implant nests
    inside the previous one (IMP.E.1-4), contacts land on a contact implant (V1.I.2) and M1
    encloses every via by more than 100 nm (M1.E.1).
    """
    c.add("SLAB", pf.Rectangle((0, -slab_half), (length, slab_half)))
    c.add("SI", pf.Rectangle((0, -core_half), (length, core_half)))
    c.add("PW", pf.Rectangle((0, -slab_half), (length, 0.06)))
    c.add("NW", pf.Rectangle((0, 0.06), (length, slab_half)))
    c.add("PP", pf.Rectangle((0, -slab_half + 0.05), (length, -core_half)))
    c.add("NP", pf.Rectangle((0, core_half), (length, slab_half - 0.05)))
    c.add("PPP", pf.Rectangle((0, -slab_half + 0.10), (length, -slab_half + 1.05)))
    c.add("NPP", pf.Rectangle((0, slab_half - 1.05), (length, slab_half - 0.10)))
    y_p = -slab_half + 0.575
    y_n = slab_half - 0.575
    n_via = max(1, int((length - 1.0) // via_pitch))
    for i in range(n_via):
        x = 0.5 + i * via_pitch
        for y in (y_p, y_n):
            c.add("VIA1", pf.Rectangle(center=(x + 0.15, y), size=(0.3, 0.3)))
    for y in (y_p, y_n):
        c.add("M1", pf.Rectangle((0.2, y - 0.35), (length - 0.2, y + 0.35)))
    return y_p, y_n


@pf.parametric_component(name_prefix="pnps")
def pn_phase_shifter(*, length: pft.PositiveDimension = 200.0, technology=None):
    """Rib waveguide with a lateral depletion-mode p-n junction.

    The electro-optic strength is the MSB segment of the two-segment Z-shape junction in
    Yuan et al. 2024: V_pi.L = 0.53 V.cm.
    """
    tech = _tech(technology)
    spec = _spec(tech, RIB)
    c = pf.Component(technology=tech)
    y_p, y_n = _junction_rows(c, length)
    c.add_port(pf.Port((0, 0), 0, spec), "P0")
    c.add_port(pf.Port((length, 0), 180, spec), "P1")
    c.add_terminal({
        "ANODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((0.2, y_p - 0.35), (4.0, y_p + 0.35))),
        "CATHODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle((0.2, y_n - 0.35), (4.0, y_n + 0.35))),
    })
    c.add_model(_wg_model(length, port_spec=RIB, v_piL=cal.V_PI_L, voltage=0.0), "analytic")
    c.activate_model("analytic")
    _named(c, "pnps", l=length)
    return add_recognition(c, icon="eo_ps")


@pf.parametric_component(name_prefix="mrm")
def micro_ring_modulator(*, ring_radius: pft.PositiveDimension = 11.9752,
                         coupling_gap: pft.PositiveDimension = 0.18,
                         coupling_length: pft.PositiveDimension = 4.0,
                         euler_fraction: float = 0.15,
                         s_bend_offset: pft.PositiveDimension = 5.0,
                         taper_length: pft.PositiveDimension = 10.0,
                         heater_width: pft.PositiveDimension = 1.2,
                         voltage: pft.Voltage = 0.0,
                         draw_junction: bool = True, technology=None):
    """Depletion-mode microring modulator with a radial p-n junction and a TiN heater.

    The ring is a true circle, so its optical round trip is exactly `2 pi R` and the model
    reads its length from the geometry. That matters for a WDM bank: the waveguide index is
    a fixed property of the cross-section, and it is the radius that sets which resonance a
    channel lands on. The default radius puts mode 160 at 1310.0 nm.

    `coupling_length` is the straight run of the bus alongside the ring; it sets the coupling
    strength, not the ring circumference.

    `voltage` is a cell parameter rather than something you set on the model afterwards: the
    S-matrix is cached per component, so mutating a model in place returns the stale result.

    The junction, contacts and heater follow the device in Yuan et al. 2024; the model is a
    RingModel paired with a RingTimeStepper so the same cell drives a time-domain link.
    """
    tech = _tech(technology)
    spec = _spec(tech, RIB)
    w, _ = _profile(spec, tech, "SI")
    c = pf.Component(technology=tech)
    s_bend_length = 5 * s_bend_offset

    def ring_band(width):
        """Annular band of the given width, centred on the ring centreline."""
        return pf.Circle(radius=ring_radius + width / 2, inner_radius=ring_radius - width / 2)

    ring = ring_band(w)
    y_bus = -(ring_radius + w + coupling_gap)  # bus centreline, tangent below the ring
    x_port = -coupling_length / 2 - s_bend_length - taper_length
    y_port = y_bus - s_bend_offset
    bus = (pf.Path((x_port, y_port), w)
           .segment((taper_length, 0), relative=True)
           .s_bend((s_bend_length, s_bend_offset), euler_fraction, relative=True)
           .segment((coupling_length, 0), relative=True)
           .s_bend((s_bend_length, -s_bend_offset), euler_fraction, relative=True)
           .segment((-x_port, y_port)))
    c.add("SI", ring, bus)

    filled = pf.envelope(ring, offset=0.0)  # solid disc out to the core's outer edge

    def reg(d):
        """Filled region whose boundary sits `d` from the ring centreline (d may be < 0)."""
        out = pf.envelope(ring, offset=d - 0.25) if d >= 0.25 else pf.offset(filled, d - 0.25)
        return out if isinstance(out, list) else [out]

    def opened(shapes, r=0.2):
        """Morphological opening: drop fragments thinner than 2r left by the keep-out cut.

        `round_joins` matters: with the default miter joins the outward pass rebuilds sharp
        corners and leaves wedges narrower than the metal width rule.
        """
        return pf.offset(pf.offset(shapes, -r, round_joins=True), r, round_joins=True)

    def closed(shapes, r=0.3):
        """Morphological closing: merge regions separated by less than 2r."""
        return pf.offset(pf.offset(shapes, r, round_joins=True), -r, round_joins=True)

    bus_keepout = pf.envelope(bus, offset=0.6)  # the coupling region stays undoped

    def outer(d_hi, d_lo):
        return opened(pf.boolean(pf.boolean(reg(d_hi), reg(d_lo), "-"), bus_keepout, "-"))

    # trim the envelope at the x extremes, or the slab runs 1.75 um past the cell's own
    # optical ports and every butt-connected transition overlaps this device
    bus_slab = pf.envelope(bus, offset=1.75, trim_x_min=True, trim_x_max=True)
    slab = closed(pf.boolean([*reg(2.0), bus_slab], [], "+"), 0.6)
    c.add("SLAB", *slab)

    if draw_junction:
        # Nesting is enforced by construction: each graded implant is intersected with the
        # one below it, so IMP.E.1-E.4 cannot be broken by the keep-out cut or by opening.
        p_well = opened(pf.boolean(reg(-0.06), bus_keepout, "-"))
        p_imp = pf.boolean(opened(pf.boolean(reg(-0.25), bus_keepout, "-")), p_well, "*")
        p_con = pf.boolean(opened(pf.boolean(reg(-0.70), bus_keepout, "-")), p_imp, "*")
        n_well = outer(1.75, 0.06)
        n_imp = pf.boolean(outer(1.70, 0.25), n_well, "*")
        n_con = pf.boolean(outer(1.65, 0.70), n_imp, "*")
        c.add("PW", *p_well, "NW", *n_well, "PP", *p_imp, "NP", *n_imp,
              "PPP", *p_con, "NPP", *n_con)
        # 1.6 um wide so a via landing fits inside the band without a notch
        m1_in = opened(pf.boolean(pf.boolean(reg(-0.60), reg(-2.20), "-"), bus_keepout, "-"), 0.3)
        m1_out = opened(pf.boolean(pf.boolean(reg(1.60), reg(0.80), "-"), bus_keepout, "-"), 0.3)
        c.add("M1", *m1_in, *m1_out)

        # Contact rings. A via is kept only if it is fully covered by both the contact
        # implant and M1, which is exactly what V1.I.2 and M1.E.1 ask for; near the coupling
        # region those layers have been cut away, so those vias simply do not get placed.
        def covered_by(rect, shapes, margin=0.0):
            probe = pf.Rectangle(center=rect.center, size=(rect.size[0] + 2 * margin,
                                                           rect.size[1] + 2 * margin))
            return len(pf.boolean([probe], shapes, "-")) == 0

        n_via = 44
        for i in range(n_via):
            a = 2 * np.pi * i / n_via
            # the two heater leads cross the outer contact ring on their way out, and VIA1
            # reaches through the heater level, so the vias they would land on are omitted
            under_lead = min(abs((np.degrees(a) - t + 180) % 360 - 180) for t in _LEAD_ANGLES)
            for sign, con, m1 in ((-1, p_con, m1_in), (1, n_con, m1_out)):
                if sign > 0 and under_lead < 6.0:
                    continue
                r = ring_radius + sign * 1.2
                rect = pf.Rectangle(center=(r * np.cos(a), r * np.sin(a)), size=(0.3, 0.3))
                if covered_by(rect, con) and covered_by(rect, m1, margin=0.12):
                    c.add("VIA1", rect)
        c.add_terminal({
            "ANODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle(center=(0, ring_radius - 1.2), size=(2.0, 0.5))),
            "CATHODE": pf.Terminal(_layer(tech, "M1"), pf.Rectangle(center=(0, ring_radius + 1.2), size=(2.0, 0.5))),
        })

    # TiN heater: an arc that stops just short of a full turn, with one radial lead at each
    # end. Three things drive the shape. It has to be broken, because a closed ring with two
    # taps just shorts the short way round and heats a quarter of the device. The leads come
    # out at 45 and 135 degrees rather than straight up, which leaves the top of the ring
    # clear for the cathode to escape northwards to its own pad. And the band is narrower
    # than the contact pitch so the heater keeps its distance from both via rings: VIA1 runs
    # from the silicon all the way up to M1, straight through the heater level, so heater
    # metal overlapping a contact via is a short, not a spacing nicety.
    #
    # Lead, arc and lead are ONE path rather than a union of three shapes. Unioning a
    # straight lead into a polygon-approximated arc leaves sub-100 nm nicks and spikes where
    # the two boundaries, built from different chord decompositions, cross at a shallow
    # angle, and cleaning those up morphologically just moves the artifacts onto the arc's
    # own chord vertices. A single path has the joins computed once, exactly.
    def _u(angle):
        return np.array([np.cos(np.radians(angle)), np.sin(np.radians(angle))])

    a_plus, a_minus = _LEAD_ANGLES
    lead_end = ring_radius + 7.2
    # The arc and the two leads are separate shapes that overlap, not one path with joins in
    # it. A single path reads better and asks the path builder to intersect a straight section
    # with a curved one at each lead, which below about 1 nm of curve tolerance it cannot:
    # "Unable to find intersection between path sections", for one radius out of five. That
    # is the kind of failure that arrives the day somebody sets a finer tolerance for an
    # unrelated reason - the write-back flow does exactly that - so there is no join to solve.
    #
    # What overlapping shapes cost is a tangential wedge where each lead meets the arc, whose
    # tip is narrower than any space rule. The deck answers that with angle_limit(80) on the
    # two HTR rules rather than by deforming the geometry; see the comment there.
    c.add("HTR", pf.Circle(radius=ring_radius + heater_width / 2,
                           inner_radius=ring_radius - heater_width / 2,
                           sector=(a_plus, a_minus + 360)))
    for angle in _LEAD_ANGLES:
        u = _u(angle)
        n = np.array([-u[1], u[0]]) * (heater_width / 2)
        inner = u * (ring_radius - heater_width)
        outer = u * lead_end
        c.add("HTR", pf.Polygon([tuple(float(v) for v in inner - n),
                                 tuple(float(v) for v in outer - n),
                                 tuple(float(v) for v in outer + n),
                                 tuple(float(v) for v in inner + n)]))

    # HTRP/HTRN, not H+/H-. A terminal name ends up as a net name in the SPICE netlist a
    # write-back emits, and spectre reads `+` as an operator: the netlist parses as far as
    # `e_I2_H+` and stops, spiceIn creates no schematic, and Check Against Source has nothing
    # to check. Cosmetic-looking names are part of the interface.
    leads = {}
    for name, angle in zip(("HTRP", "HTRN"), _LEAD_ANGLES):
        # Snapped to the layout grid, and deliberately wider than the 4 um trace that lands
        # on it. Both matter. A pad whose edge falls half a grid step from the route's edge
        # rounds the other way and leaves a 1 nm sliver that the width rule finds, and a pad
        # exactly as wide as its route has two edges that can disagree that way.
        tip = tuple(round(float(v), 3) for v in _u(angle) * (ring_radius + 6.0))
        c.add("VIA2", pf.Rectangle(center=tip, size=(0.4, 0.4)))
        c.add("M1", pf.Rectangle(center=tip, size=(5.2, 5.2)))
        leads[name] = pf.Terminal(_layer(tech, "M1"), pf.Rectangle(center=tip, size=(5.2, 5.2)))
    c.add_terminal(leads)

    c.add_port(pf.Port((x_port, y_port), 0, spec), "P0")
    c.add_port(pf.Port((-x_port, y_port), 180, spec), "P1")
    # The junction drive, as a port rather than only as the ANODE/CATHODE terminals. The
    # terminals are the contacts a router lands metal on; this is where a circuit simulator
    # applies a voltage, and `RingModel`/`RingTimeStepper` read it to move the resonance. A
    # modulator with no drive port can be laid out and verified and cannot be simulated, which
    # is half a model. It draws nothing and takes no pin - see `_recognition.has_pin_marker`.
    c.add_port(pf.Port((0.0, ring_radius + 6.0), -90,
                       pf.virtual_port_spec(classification="electrical")), "DRIVE")

    # The optical round trip comes from the geometry, and the index is a property of the
    # cross-section, not something re-tuned per cell: that is what makes the radius the knob
    # that sets a WDM channel.
    length = 2 * np.pi * ring_radius
    model = pf.RingModel(kappa1=cal.KAPPA1, kappa2=None, n_eff=cal.N_EFF_RIB,
                         n_group=cal.N_GROUP_RIB, length=length,
                         propagation_loss=cal.ALPHA0_DBUM,
                         reference_frequency=cal.REF_FREQUENCY, voltage=voltage,
                         dn_dv=cal.DN_DV)
    model.time_stepper = pf.RingTimeStepper(
        kappa1=cal.KAPPA1, kappa2=None, n_eff=cal.N_EFF_RIB, n_group=cal.N_GROUP_RIB,
        length=length, propagation_loss=cal.ALPHA0_DBUM, f_3dB=cal.F_3DB_RC,
        dn_dv=cal.DN_DV)
    c.add_model(model, "ring")
    c.activate_model("ring")
    _named(c, "mrm", r=ring_radius, v=voltage)
    return add_recognition(c, icon="mrm")


# =====================================================================
# Tier 2: the nitride layer, and passives the transmitter does not need
# =====================================================================
# A SiN waveguide and a SiN bend are not new cells: `straight(port_spec=SIN)` and
# `bend(port_spec=SIN, radius=...)` already draw on the nitride layer and, since `_GUIDE`
# exists, are modelled with the nitride index. What nitride did need was a way on and off it,
# which is the two cells below.


@pf.parametric_component(name_prefix="s2n")
def si_to_sin(*, length: pft.PositiveDimension = 40.0,
              si_tip: pft.PositiveDimension = 0.15,
              sin_tip: pft.PositiveDimension = 0.30, technology=None):
    """Interlayer transition: silicon down to a tip, nitride up from a tip, overlapping.

    Both guides taper across the same length in opposite directions, so the mode is handed
    over adiabatically rather than launched across a gap. `length` is the whole overlap; the
    default is long because the two indices are far apart - 2.59 against 1.64 - and a short
    transition radiates rather than converts.

    Both tips sit exactly on the deck's minimum width - 150 nm on SI.W.1, 300 nm on SIN.W.1 -
    rather than below it. A narrower tip converts better and this PDK's rule is to design to
    its own deck rather than to waive it, so the tip is as sharp as the process allows and no
    sharper. Drawn at 250 nm, which is what an inverse taper would like, this cell reports
    SIN.W.1 and is correct to.
    """
    tech = _tech(technology)
    si_spec, sin_spec = _spec(tech, STRIP), _spec(tech, SIN)
    si_w, _ = _profile(si_spec, tech, "SI")
    sin_w, _ = _profile(sin_spec, tech, "SIN")
    c = pf.Component(technology=tech)
    c.add("SI", pf.Polygon([(0, si_w / 2), (length, si_tip / 2),
                            (length, -si_tip / 2), (0, -si_w / 2)]))
    c.add("SIN", pf.Polygon([(0, sin_tip / 2), (length, sin_w / 2),
                             (length, -sin_w / 2), (0, -sin_tip / 2)]))
    c.add_port(pf.Port((0, 0), 0, si_spec), "P0")
    c.add_port(pf.Port((length, 0), 180, sin_spec), "P1")

    # Half the length in each guide, which is what an adiabatic hand-over is to first order.
    # The mean group index is the honest single number for a cell whose index changes along
    # it; anything needing better than that wants EME, not a compact model.
    n_eff = 0.5 * (cal.N_EFF_STRIP + cal.N_EFF_SIN)
    n_g = 0.5 * (cal.N_GROUP_STRIP + cal.N_GROUP_SIN)
    c.add_model(pf.AnalyticWaveguideModel(
        n_eff=n_eff, n_group=n_g, length=length,
        propagation_loss=0.5 * (cal.LOSS_STRIP + cal.LOSS_SIN),
        extra_loss=cal.LOSS_SI_TO_SIN, reference_frequency=cal.REF_FREQUENCY), "analytic")
    c.activate_model("analytic")
    _named(c, "s2n", l=length)
    return add_recognition(c, icon="transition")


@pf.parametric_component(name_prefix="ec")
def sin_edge_coupler(*, length: pft.PositiveDimension = 150.0,
                     tip_width: pft.PositiveDimension = 0.30,
                     setback: pft.PositiveDimension = 5.0, technology=None):
    """Inverse-taper edge coupler on nitride, for a lensed fibre at the die facet.

    Nitride rather than silicon because the point of an edge coupler is to meet a fibre mode,
    and a 1 um nitride guide is already most of the way there: the taper only has to expand
    the mode the rest of the way, instead of fighting silicon's confinement.

    The facet port is a `FiberPort` like the grating coupler's, but pointing sideways rather
    than up - the fibre arrives in the plane of the die. `setback` is the gap between the tip
    and the port, which is where the dicing street goes; it is drawn as nothing, because the
    edge of the chip is an absence.

    The tip is 300 nm because that is SIN.W.1's minimum. A real nitride edge coupler would
    want it narrower and would carry a process waiver to get there; this PDK designs to its
    deck instead, and the cost is a fraction of a dB of mode mismatch.
    """
    tech = _tech(technology)
    spec = _spec(tech, SIN)
    w, _ = _profile(spec, tech, "SIN")
    c = pf.Component(technology=tech)
    # Tip at x = 0 facing the facet, full width at x = length facing the circuit.
    c.add("SIN", pf.Polygon([(0, tip_width / 2), (length, w / 2),
                             (length, -w / 2), (0, -tip_width / 2)]))
    c.add_port(pf.Port((length, 0), 180, spec), "P0")
    # Direction (-1, 0, 0): straight out of the facet, in plane, no tilt. The mode field
    # diameter is the lensed-fibre spot this taper is designed to meet.
    c.add_port(pf.FiberPort((-setback, 0, 0), (-1, 0, 0), (3.2, 3.2), (-2.0, 2.0)), "FACET")
    c.add_model(_edge_coupler_model(("P0", "FACET")), "measured")
    c.activate_model("measured")
    _named(c, "ec", l=length)
    return add_recognition(c, icon="edge_coupler")


def _edge_coupler_model(port_names):
    """Synthetic measurement data for the edge coupler.

    Broadband and gently sloped, which is how an inverse taper behaves and the reason a PDK
    offers one next to a grating coupler: the grating is 20 nm wide and this is flat across
    the band. Same `DataModel` path a measured .s2p would take.
    """
    lam = np.linspace(1.26, 1.36, 101)
    loss_db = cal.EC_PEAK_LOSS_DB + cal.EC_SLOPE_DB_PER_NM * (lam - cal.WAVELENGTH) * 1e3
    t = 10 ** (-loss_db / 20)
    r = 10 ** (-30 / 20)
    freqs = pf.C_0 / lam
    order = np.argsort(freqs)
    s_array = np.zeros((lam.size, 2, 2), dtype=complex)
    s_array[:, 0, 0] = r
    s_array[:, 1, 1] = r
    s_array[:, 0, 1] = t[order]
    s_array[:, 1, 0] = t[order]
    return pf.DataModel(s_array=s_array, frequencies=freqs[order], ports=list(port_names))


@pf.parametric_component(name_prefix="mrr")
def micro_ring_resonator(*, ring_radius: pft.PositiveDimension = 11.9942,
                         coupling_gap: pft.PositiveDimension = 0.18,
                         coupling_length: pft.PositiveDimension = 4.0,
                         heater_width: pft.PositiveDimension = 1.2, technology=None):
    """Passive ring resonator with a heater, and no junction.

    The same racetrack as the modulator with the doping left out, which is what a filter or a
    sensor wants: no junction means no free-carrier loss and a higher Q, and the heater is
    still there because a passive ring that cannot be trimmed onto its channel is not usable
    on a real die.
    """
    c = micro_ring_modulator(ring_radius=ring_radius, coupling_gap=coupling_gap,
                             coupling_length=coupling_length, heater_width=heater_width,
                             draw_junction=False, technology=technology)
    # The modulator's drive port goes with its junction. Leaving it on would advertise an
    # electrical input that moves nothing, and a circuit built against it would apply a
    # voltage and see no phase.
    c.remove_port("DRIVE")

    # A fresh model rather than the modulator's: with no junction there is no dn_dv, and
    # leaving it in would let a voltage tune a ring that has nothing to tune it with.
    length = 2 * np.pi * ring_radius
    model = pf.RingModel(kappa1=cal.KAPPA1, kappa2=None, n_eff=cal.N_EFF_RIB,
                         n_group=cal.N_GROUP_RIB, length=length,
                         propagation_loss=cal.LOSS_RIB,
                         reference_frequency=cal.REF_FREQUENCY,
                         dn_dT=cal.DN_DT, temperature=293.0)
    model.time_stepper = pf.RingTimeStepper(
        kappa1=cal.KAPPA1, kappa2=None, n_eff=cal.N_EFF_RIB, n_group=cal.N_GROUP_RIB,
        length=length, propagation_loss=cal.LOSS_RIB, dn_dT=cal.DN_DT)
    for name in list(c.models):
        c.remove_model(name)
    c.add_model(model, "ring")
    c.activate_model("ring")
    _named(c, "mrr", r=ring_radius)
    # NOT `add_recognition` again. The modulator this was built from already drew its DEVREC
    # and its pin markers, and a second DEVREC on top of the first is two overlapping device
    # boundaries, which DEV.1 reports and is right to. Only the palette icon changes.
    c.properties["__thumbnail__"] = "mrr"
    return c
