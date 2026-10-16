"""SiPhoPDK layer table.

SiPhoPDK is a SYNTHETIC demonstration PDK. Layer numbers, names and purposes are our own
invention. They are not taken from any foundry PDK.
The stack it represents (220 nm SOI + partial-etch slab + SiN + Ge + heater + 2 metals) is
the generic O-band datacom platform described in the public literature, e.g.
S. K. Yeh et al., "Silicon Photonics Platform for Next Generation Data Communication
Technologies", Fig. 1b and Fig. 5.
"""

import photonforge as pf

# min_width and min_spacing repeat the width and space rules of rules/siphopdk.drc, so that
# PhotonForge routers draw to them. `conductive` marks the layers that carry current.

_layers = {
    # ---- optical structural layers -------------------------------------------------
    "SI": pf.LayerSpec((1, 0), "Silicon full etch, 220 nm core", "#78838fc0", "solid",
                       min_width=0.15, min_spacing=0.15),
    "SLAB": pf.LayerSpec((2, 0), "Silicon partial etch, 90 nm slab remains", "#98a3ae90", "//",
                         min_width=0.30, min_spacing=0.30),
    "GRAT": pf.LayerSpec((3, 0), "Grating shallow etch, 70 nm", "#ffd16690", "\\\\",
                         min_width=0.10, min_spacing=0.10),
    "SIN": pf.LayerSpec((5, 0), "Silicon nitride, 300 nm", "#2f5fd0c0", "solid",
                        min_width=0.30, min_spacing=0.30),
    "GE": pf.LayerSpec((8, 0), "Germanium epi, 250 nm", "#38b0d4c0", "solid",
                       min_width=0.60, min_spacing=1.00),
    # ---- implants ------------------------------------------------------------------
    "NW": pf.LayerSpec((20, 0), "n well, 3e18 cm-3", "#bfe3b570", ".", min_width=0.30),
    "PW": pf.LayerSpec((21, 0), "p well, 3e18 cm-3", "#f3b6c870", ".", min_width=0.30),
    "NP": pf.LayerSpec((22, 0), "n+ implant, 6e18 cm-3", "#8ecb7e70", ":"),
    "PP": pf.LayerSpec((23, 0), "p+ implant, 6e18 cm-3", "#e8809f70", ":"),
    "NPP": pf.LayerSpec((24, 0), "n++ contact implant, 2e20 cm-3", "#4f9e4290", "xx",
                        conductive=True),
    "PPP": pf.LayerSpec((25, 0), "p++ contact implant, 9e19 cm-3", "#c94f7490", "xx",
                        conductive=True),
    # ---- heater and metals ---------------------------------------------------------
    # The heater is a resistor, not a wire, so it is not conductive: its two leads are two
    # nets. Silicon is not conductive either; only the contact implants NPP and PPP conduct.
    "HTR": pf.LayerSpec((30, 0), "TiN heater, 120 nm", "#c0272dc0", "x",
                        min_width=0.50, min_spacing=1.00),
    "VIA1": pf.LayerSpec((40, 0), "Contact via, silicon/Ge to M1", "#8a5a2bc0", "solid",
                         conductive=True, min_width=0.18, min_spacing=0.24),
    "VIA2": pf.LayerSpec((41, 0), "Via, heater to M1", "#8a5a2bc0", "solid", conductive=True),
    "VIA3": pf.LayerSpec((42, 0), "Via, M1 to M2", "#a3652fc0", "solid",
                         conductive=True, min_width=0.25, min_spacing=0.35),
    "M1": pf.LayerSpec((50, 0), "Metal 1, 500 nm", "#c9793fb0", "\\\\",
                       conductive=True, min_width=0.50, min_spacing=0.50),
    "M2": pf.LayerSpec((52, 0), "Metal 2, 1000 nm, RF and power", "#f0a71fb0", "//",
                       conductive=True, min_width=1.00, min_spacing=1.00),
    "PAD": pf.LayerSpec((60, 0), "Passivation opening, bond pad", "#ffd166b0", "solid",
                        min_width=40.0),
    # ---- process and mask prep -----------------------------------------------------
    "DT": pf.LayerSpec((70, 0), "Deep trench, edge coupler facet", "#7f9db960", "hollow"),
    "FLOORPLAN": pf.LayerSpec((80, 0), "Usable floorplan area", "#00643C50", "hollow"),
    "DICE": pf.LayerSpec((81, 0), "Dicing lane", "#33333350", "hollow"),
    "FILL_SI": pf.LayerSpec((90, 0), "Silicon density fill", "#78838f40", "."),
    "FILL_M": pf.LayerSpec((91, 0), "Metal density fill", "#f0a71f40", "."),
    "NOFILL": pf.LayerSpec((92, 0), "Fill exclusion", "#88888840", "hollow"),
    "EXCL_DRC": pf.LayerSpec((95, 0), "DRC waiver region", "#ff000040", "hollow"),
    # ---- recognition layers: the SiPhoPDK LVS convention ---------------------------
    # Every PDK cell draws exactly one DEVREC boundary, one PINREC marker per optical
    # port and one PINREC_E marker per electrical terminal, each with a TEXT label
    # carrying the port name. Extraction reads devices from DEVREC, nets from the
    # connectivity table, and port identity from PINREC/PINREC_E + TEXT.
    "DEVREC": pf.LayerSpec((100, 0), "Device recognition boundary", "#00963050", "hollow"),
    "PINREC": pf.LayerSpec((101, 0), "Optical pin marker", "#00c85380", "solid"),
    "PINREC_E": pf.LayerSpec((102, 0), "Electrical pin marker", "#ff8c0080", "solid"),
    "TEXT": pf.LayerSpec((105, 0), "Labels, not fabricated", "#20202080", "hollow"),
    "NETLBL": pf.LayerSpec((106, 0), "Net label, not fabricated", "#3355ff80", "hollow"),
}
