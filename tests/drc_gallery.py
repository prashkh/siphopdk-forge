"""The SiPhoPDK DRC violation gallery: one deliberate violation per rule.

Each rule id gets its own site containing the smallest geometry that breaks exactly that
rule (helper geometry is drawn generously so a site does not trip unrelated rules). The
gallery is the acceptance test for the deck: every rule must fire on its own site, and the
companion `clean` cell must produce no violations at all.
"""

import photonforge as pf

PITCH = 70.0
COLS = 7


def R(c, layer, x0, y0, x1, y1, dx=0.0, dy=0.0):
    c.add(layer, pf.Rectangle((x0 + dx, y0 + dy), (x1 + dx, y1 + dy)))


# Each builder draws its violation in a local frame around (0, 0).
BUILDERS = {}


def rule(rid):
    def deco(fn):
        BUILDERS[rid] = fn
        return fn

    return deco


# ---- A. flow ---------------------------------------------------------------------
@rule("FP.1")
def _(c, dx, dy):
    R(c, "SI", 0, 0, 2, 2, dx, dy)  # placed outside the floorplan by the caller


@rule("DEV.1")
def _(c, dx, dy):
    R(c, "DEVREC", 0, 0, 6, 6, dx, dy)
    R(c, "DEVREC", 4, 4, 10, 10, dx, dy)


# ---- B. silicon ------------------------------------------------------------------
@rule("SI.W.1")
def _(c, dx, dy):
    R(c, "SI", 0, 0, 5, 0.10, dx, dy)


@rule("SI.S.1")
def _(c, dx, dy):
    R(c, "SI", 0, 0, 5, 1, dx, dy)
    R(c, "SI", 0, 1.10, 5, 2.10, dx, dy)


@rule("SI.A.1")
def _(c, dx, dy):
    R(c, "SI", 0, 0, 0.17, 0.17, dx, dy)  # legal width, area 0.029 um^2


# ---- C. slab ---------------------------------------------------------------------
@rule("SLAB.W.1")
def _(c, dx, dy):
    R(c, "SLAB", 0, 0, 5, 0.20, dx, dy)


@rule("SLAB.S.1")
def _(c, dx, dy):
    R(c, "SLAB", 0, 0, 5, 2, dx, dy)
    R(c, "SLAB", 0, 2.20, 5, 4.20, dx, dy)


@rule("SLAB.E.1")
def _(c, dx, dy):
    R(c, "SLAB", 0, 0, 6, 0.9, dx, dy)
    R(c, "SI", 0, 0.2, 6, 0.7, dx, dy)  # only 200 nm of slab on each side


# ---- D. grating ------------------------------------------------------------------
@rule("GRAT.W.1")
def _(c, dx, dy):
    R(c, "SI", 0, 0, 6, 6, dx, dy)
    R(c, "GRAT", 1, 1, 5, 1.05, dx, dy)


@rule("GRAT.S.1")
def _(c, dx, dy):
    R(c, "SI", 0, 0, 6, 6, dx, dy)
    R(c, "GRAT", 1, 1, 5, 2, dx, dy)
    R(c, "GRAT", 1, 2.05, 5, 3, dx, dy)


@rule("GRAT.I.1")
def _(c, dx, dy):
    R(c, "GRAT", 0, 0, 3, 3, dx, dy)  # no silicon underneath


# ---- E. nitride ------------------------------------------------------------------
@rule("SIN.W.1")
def _(c, dx, dy):
    R(c, "SIN", 0, 0, 5, 0.20, dx, dy)


@rule("SIN.S.1")
def _(c, dx, dy):
    R(c, "SIN", 0, 0, 5, 2, dx, dy)
    R(c, "SIN", 0, 2.20, 5, 4.20, dx, dy)


# ---- F. germanium ----------------------------------------------------------------
@rule("GE.W.1")
def _(c, dx, dy):
    R(c, "SLAB", 0, 0, 6, 4, dx, dy)
    R(c, "GE", 1, 1, 5, 1.40, dx, dy)


@rule("GE.S.1")
def _(c, dx, dy):
    R(c, "SLAB", 0, 0, 10, 8, dx, dy)
    R(c, "GE", 1, 1, 4, 4, dx, dy)
    R(c, "GE", 4.5, 1, 8, 4, dx, dy)


@rule("GE.E.1")
def _(c, dx, dy):
    R(c, "SLAB", 0, 0, 6, 6, dx, dy)
    R(c, "GE", 0.1, 0.1, 5.9, 5.9, dx, dy)  # only 100 nm of slab around it


# ---- G. implants -----------------------------------------------------------------
@rule("IMP.W.1")
def _(c, dx, dy):
    R(c, "NW", 0, 0, 5, 0.20, dx, dy)


@rule("IMP.O.1")
def _(c, dx, dy):
    R(c, "NW", 0, 0, 4, 4, dx, dy)
    R(c, "PW", 2, 2, 6, 6, dx, dy)


@rule("IMP.E.1")
def _(c, dx, dy):
    R(c, "NP", 0, 0, 3, 3, dx, dy)  # no NW under it


@rule("IMP.E.2")
def _(c, dx, dy):
    R(c, "PP", 0, 0, 3, 3, dx, dy)


@rule("IMP.E.3")
def _(c, dx, dy):
    R(c, "SLAB", -1, -1, 4, 4, dx, dy)
    R(c, "NPP", 0, 0, 3, 3, dx, dy)  # on silicon, but no NP around it


@rule("IMP.E.4")
def _(c, dx, dy):
    R(c, "SLAB", -1, -1, 4, 4, dx, dy)
    R(c, "PPP", 0, 0, 3, 3, dx, dy)


@rule("IMP.I.1")
def _(c, dx, dy):
    R(c, "NW", -1, -1, 4, 4, dx, dy)
    R(c, "NP", -0.5, -0.5, 3.5, 3.5, dx, dy)
    R(c, "NPP", 0, 0, 3, 3, dx, dy)  # nested correctly, but no silicon under it


# ---- H. contacts, heater, metals --------------------------------------------------
def _good_contact(c, dx, dy, x0, y0, size=0.25):
    """A compliant VIA1 stack: slab + NPP + NP + NW + via + enclosing M1."""
    R(c, "SLAB", x0 - 1, y0 - 1, x0 + size + 1, y0 + size + 1, dx, dy)
    R(c, "NW", x0 - 0.9, y0 - 0.9, x0 + size + 0.9, y0 + size + 0.9, dx, dy)
    R(c, "NP", x0 - 0.8, y0 - 0.8, x0 + size + 0.8, y0 + size + 0.8, dx, dy)
    R(c, "NPP", x0 - 0.7, y0 - 0.7, x0 + size + 0.7, y0 + size + 0.7, dx, dy)


@rule("V1.W.1")
def _(c, dx, dy):
    _good_contact(c, dx, dy, 0, 0)
    R(c, "VIA1", 0, 0, 0.15, 0.15, dx, dy)
    R(c, "M1", -0.6, -0.6, 0.75, 0.75, dx, dy)


@rule("V1.S.1")
def _(c, dx, dy):
    _good_contact(c, dx, dy, 0, 0, size=1.5)
    R(c, "VIA1", 0, 0, 0.25, 0.25, dx, dy)
    R(c, "VIA1", 0.45, 0, 0.70, 0.25, dx, dy)  # 200 nm apart
    R(c, "M1", -0.6, -0.6, 1.3, 0.85, dx, dy)


@rule("V1.I.1")
def _(c, dx, dy):
    R(c, "VIA1", 0, 0, 0.3, 0.3, dx, dy)  # lands on nothing
    R(c, "M1", -0.6, -0.6, 0.9, 0.9, dx, dy)


@rule("V1.I.2")
def _(c, dx, dy):
    R(c, "SLAB", -1, -1, 1.3, 1.3, dx, dy)
    R(c, "VIA1", 0, 0, 0.3, 0.3, dx, dy)  # on silicon, no contact implant
    R(c, "M1", -0.6, -0.6, 0.9, 0.9, dx, dy)


@rule("M1.W.1")
def _(c, dx, dy):
    R(c, "M1", 0, 0, 6, 0.30, dx, dy)


@rule("M1.S.1")
def _(c, dx, dy):
    R(c, "M1", 0, 0, 6, 2, dx, dy)
    R(c, "M1", 0, 2.30, 6, 4.30, dx, dy)


@rule("M1.E.1")
def _(c, dx, dy):
    _good_contact(c, dx, dy, 0, 0)
    R(c, "VIA1", 0, 0, 0.25, 0.25, dx, dy)
    R(c, "M1", -0.05, -0.05, 2.5, 2.5, dx, dy)  # wide metal, but only 50 nm past two via edges


@rule("M1.E.2")
def _(c, dx, dy):
    R(c, "HTR", -1, -1, 1.3, 1.3, dx, dy)
    R(c, "VIA2", 0, 0, 0.3, 0.3, dx, dy)
    R(c, "M1", -0.05, -0.05, 2.5, 2.5, dx, dy)


@rule("HTR.W.1")
def _(c, dx, dy):
    R(c, "HTR", 0, 0, 6, 0.30, dx, dy)


@rule("HTR.S.1")
def _(c, dx, dy):
    R(c, "HTR", 0, 0, 6, 2, dx, dy)
    R(c, "HTR", 0, 2.50, 6, 4.50, dx, dy)


@rule("V2.I.1")
def _(c, dx, dy):
    R(c, "VIA2", 0, 0, 0.3, 0.3, dx, dy)  # no heater under it
    R(c, "M1", -0.6, -0.6, 0.9, 0.9, dx, dy)


def _good_via3(c, dx, dy, x0, y0, size=0.3):
    R(c, "M1", x0 - 0.8, y0 - 0.8, x0 + size + 0.8, y0 + size + 0.8, dx, dy)
    R(c, "M2", x0 - 1.0, y0 - 1.0, x0 + size + 1.0, y0 + size + 1.0, dx, dy)


@rule("V3.W.1")
def _(c, dx, dy):
    _good_via3(c, dx, dy, 0, 0)
    R(c, "VIA3", 0, 0, 0.20, 0.20, dx, dy)


@rule("V3.S.1")
def _(c, dx, dy):
    _good_via3(c, dx, dy, 0, 0, size=1.0)
    R(c, "VIA3", 0, 0, 0.3, 0.3, dx, dy)
    R(c, "VIA3", 0.60, 0, 0.90, 0.3, dx, dy)  # 300 nm apart


@rule("M1.E.3")
def _(c, dx, dy):
    R(c, "M1", -0.05, -0.05, 2.5, 2.5, dx, dy)
    R(c, "M2", -1, -1, 3.5, 3.5, dx, dy)
    R(c, "VIA3", 0, 0, 0.3, 0.3, dx, dy)


@rule("M2.W.1")
def _(c, dx, dy):
    R(c, "M2", 0, 0, 8, 0.60, dx, dy)


@rule("M2.S.1")
def _(c, dx, dy):
    R(c, "M2", 0, 0, 8, 3, dx, dy)
    R(c, "M2", 0, 3.60, 8, 6.60, dx, dy)


@rule("M2.E.1")
def _(c, dx, dy):
    R(c, "M1", -1.5, -1.5, 3.5, 3.5, dx, dy)
    R(c, "M2", -0.10, -0.10, 3.0, 3.0, dx, dy)  # wide metal, but only 100 nm past two via edges
    R(c, "VIA3", 0, 0, 0.3, 0.3, dx, dy)


@rule("PAD.I.1")
def _(c, dx, dy):
    R(c, "PAD", 0, 0, 45, 45, dx, dy)  # no M2 under it


@rule("PAD.E.1")
def _(c, dx, dy):
    R(c, "M2", -1, -1, 46, 46, dx, dy)
    R(c, "PAD", 0, 0, 45, 45, dx, dy)  # only 1 um of M2 around it


@rule("PAD.W.1")
def _(c, dx, dy):
    R(c, "M2", -3, -3, 33, 33, dx, dy)
    R(c, "PAD", 0, 0, 30, 30, dx, dy)  # pad smaller than 40 um


def build_gallery(tech):
    c = pf.Component("drc_gallery", tech)
    ids = sorted(BUILDERS)
    for i, rid in enumerate(ids):
        col, row = i % COLS, i // COLS
        dx, dy = col * PITCH, -row * PITCH
        BUILDERS[rid](c, dx, dy)
        c.add(pf.Label(rid, (dx, dy - 8)))
    rows = (len(ids) + COLS - 1) // COLS
    # floorplan covers the grid; the FP.1 cell sits deliberately outside it
    c.add("FLOORPLAN", pf.Rectangle((-20, -rows * PITCH - 20), (COLS * PITCH + 20, 30)))
    fp1 = ids.index("FP.1")
    c.add("SI", pf.Rectangle((-40, 40), (-38, 42)))
    return c


def build_waiver(tech):
    """Two identical violations; the second is covered by EXCL_DRC and must be waived."""
    c = pf.Component("drc_waiver", tech)
    c.add("SI", pf.Rectangle((0, 0), (5, 0.10)))          # SI.W.1, reported
    c.add("SI", pf.Rectangle((0, 10), (5, 10.10)))        # SI.W.1, waived
    c.add("EXCL_DRC", pf.Rectangle((-1, 9), (6, 11)))
    c.add("FLOORPLAN", pf.Rectangle((-20, -20), (30, 30)))
    return c


def build_clean(tech):
    """A compliant cell: rib waveguide, contacts, heater, metal stack and a bond pad."""
    c = pf.Component("drc_clean", tech)
    # rib waveguide with a proper slab
    c.add("SLAB", pf.Rectangle((0, -1.5), (40, 1.5)))
    c.add("SI", pf.Rectangle((0, -0.25), (40, 0.25)))
    # junction implants, correctly nested and on silicon
    c.add("PW", pf.Rectangle((0, -1.4), (40, 0.06)))
    c.add("NW", pf.Rectangle((0, 0.06), (40, 1.4)))
    c.add("PP", pf.Rectangle((0, -1.3), (40, -0.25)))
    c.add("NP", pf.Rectangle((0, 0.25), (40, 1.3)))
    c.add("PPP", pf.Rectangle((1, -1.2), (39, -0.6)))
    c.add("NPP", pf.Rectangle((1, 0.6), (39, 1.2)))
    # contacts on the contact implants, enclosed by M1
    for i in range(6):
        x = 3 + i * 6
        c.add("VIA1", pf.Rectangle((x, -1.05), (x + 0.3, -0.75)))
        c.add("VIA1", pf.Rectangle((x, 0.75), (x + 0.3, 1.05)))
    c.add("M1", pf.Rectangle((2, -1.3), (38, -0.5)))
    c.add("M1", pf.Rectangle((2, 0.5), (38, 1.3)))
    # heater above, with its own via and metal
    c.add("HTR", pf.Rectangle((0, 3), (40, 4)))
    c.add("VIA2", pf.Rectangle((1, 3.3), (1.4, 3.7)))
    c.add("M1", pf.Rectangle((0.5, 2.8), (2, 4.2)))
    # M1 to M2 and a bond pad
    c.add("VIA3", pf.Rectangle((50, 0), (50.4, 0.4)))
    c.add("M1", pf.Rectangle((49, -1), (51.5, 1.5)))
    c.add("M2", pf.Rectangle((45, -5), (105, 55)))
    c.add("PAD", pf.Rectangle((50, 0), (95, 45)))
    c.add("SIN", pf.Rectangle((0, 10), (40, 11)))
    c.add("GRAT", pf.Rectangle((5, -0.2), (5.3, 0.2)))
    c.add("DEVREC", pf.Rectangle((-1, -6), (41, 12)))
    c.add("FLOORPLAN", pf.Rectangle((-20, -20), (120, 70)))
    return c
