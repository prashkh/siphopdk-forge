"""The SiPhoPDK recognition-layer convention.

Every PDK cell draws exactly one DEVREC boundary, one PINREC marker per optical port and
one PINREC_E marker per electrical terminal, each with a TEXT label carrying the port name.
Extraction reads devices from DEVREC, port identity from the pin markers plus TEXT, and
nets from the layer connectivity table. Keeping this in one helper means every cell in the
library is extractable by construction instead of by convention.
"""

import photonforge as pf

PIN_SIZE = 0.2  # um, square marker centred on the port

# The PhotonForge GUI renders `__thumbnail__` as both the library card and the schematic
# symbol. Only these names resolve to a built-in symbol; anything else silently falls back to
# the generic "Black Box" tile, so the set is validated here.
VALID_ICONS = frozenset("""
bend black_box bondpad circulator connection crossing cw_laser dc dual_mrr edge_coupler
electrical_absolute electrical_adder electrical_amplifier electrical_differential
electrical_integral electrical_mismatch electrical_multiplier electrical_termination eo_ps
filter grating_coupler isolator laser mmi1x2 mmi2x2 mrm mrr mzi mzm optical_amplifier
photodiode psgc psr s-bend source spiral taper termination tl to_ps transition u-bend wdm wg
y-splitter
""".split())


def has_pin_marker(port):
    """Whether a port of a PDK cell gets a PINREC marker, and so exists as a pin downstream.

    Two kinds of port do not. A `FiberPort` is a fibre landing above the die, carrying a z as
    well as an x and a y, with no in-plane geometry to mark. An electrical PORT is where the
    circuit simulator injects a current, which is not the same thing as a contact - the
    contacts are the cell's terminals, and those always get a marker. The photodiode's OUT is
    the only electrical port in this library; a pin marker there would put an extra member on
    the readout net in LVS.
    """
    return isinstance(port, pf.Port) and port.classification != "electrical"


def _label(c, text, position):
    c.add(pf.Label(text, position))


def add_recognition(c, *, margin=0.0, icon=None, devrec=True):
    """Add DEVREC, pin markers and labels to a finished cell, and set its thumbnail."""
    # Measure the device before adding pin markers: the markers straddle the ports, so
    # bounds taken afterwards run 100 nm past each port and every butt-connected pair then
    # overlaps by 200 nm, which DEV.1 correctly reports as overlapping devices.
    (x0, y0), (x1, y1) = c.bounds()
    for name, port in c.ports.items():
        if has_pin_marker(port):
            x, y = port.center
            c.add("PINREC", pf.Rectangle(center=(x, y), size=(PIN_SIZE, PIN_SIZE)))
            _label(c, name, (x, y))
    for name, terminal in c.terminals.items():
        # Terminal.center is a method, Port.center is a property
        x, y = terminal.center()
        c.add("PINREC_E", pf.Rectangle(center=(x, y), size=(PIN_SIZE, PIN_SIZE)))
        _label(c, name, (x, y))

    # DEVREC hugs the cell exactly. A margin looks tidier in isolation but makes every
    # butt-connected pair overlap, and DEV.1 exists precisely to catch overlapping devices.
    # Routing features such as a via stack are not devices and carry no DEVREC at all.
    if devrec:
        c.add("DEVREC", pf.Rectangle((x0 - margin, y0 - margin), (x1 + margin, y1 + margin)))

    if icon is not None:
        if icon not in VALID_ICONS:
            raise ValueError(
                f"{icon!r} is not a schematic symbol the GUI can resolve; it would render as "
                f"a Black Box tile. Pick one of: {', '.join(sorted(VALID_ICONS))}"
            )
        c.properties["__thumbnail__"] = icon
    return c
