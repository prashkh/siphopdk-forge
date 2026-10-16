"""SiPhoPDK visualisation palette.

One palette drives both views of the PDK:

* ``LayerSpec`` colours in ``_layers.py`` control the layout view (SVG, LiveViewer, GUI).
* ``td.VisualizationSpec`` attached to each medium controls the extruded cross-section
  drawn by ``pf.tidy3d_plot``.

Tidy3D colours structures by medium, not by layer, so every role that should read as a
different material gets its own medium instance (M1 and M2 are both copper, but they are
distinct objects so they plot in different colours).

Colours follow the platform cross-section convention used in the literature: silicon stays a
cool grey (core dark, slab light), nitride blue, germanium cyan, heater red, and the whole
metal stack is a warm copper-to-gold family (contacts bronze, M1 copper, VIA3 mid bronze, M2
and pads gold) so a via can never be mistaken for the silicon underneath it.
"""

import tidy3d as td

# name -> (facecolor, edgecolor, alpha)
PALETTE = {
    "substrate": ("#8d99a6", "#5b6672", 0.85),
    "oxide": ("#dfe9f3", "#b9c8d8", 0.45),
    "si_core": ("#78838f", "#3d4750", 1.0),
    "si_slab": ("#98a3ae", "#5b6672", 1.0),
    "sin": ("#2f5fd0", "#1b3a86", 1.0),
    "ge": ("#38b0d4", "#1c7794", 1.0),
    "heater": ("#c0272d", "#7c1418", 1.0),
    # metals are a warm copper-to-gold family so they never read as silicon (cool grey)
    "contact": ("#8a5a2b", "#553415", 1.0),  # VIA1/VIA2 contacts
    "m1": ("#c9793f", "#7d3f21", 1.0),
    "via": ("#a3652f", "#653a17", 1.0),  # VIA3
    "m2": ("#f0a71f", "#a86f06", 1.0),
    "opening": ("#fff3cf", "#d9b64a", 0.55),
    "trench": ("#e3eef8", "#5b7fa6", 0.65),
}

# p-type family pink, n-type family green; darker with increasing dose
IMPLANT_COLORS = {
    "PW": "#f3b6c8",
    "PP": "#e8809f",
    "PPP": "#c94f74",
    "NW": "#bfe3b5",
    "NP": "#8ecb7e",
    "NPP": "#4f9e42",
}


def viz(role):
    """VisualizationSpec for a palette role."""
    face, edge, alpha = PALETTE[role]
    return td.VisualizationSpec(facecolor=face, edgecolor=edge, alpha=alpha)


def with_viz(medium, role):
    """Attach the palette entry for ``role`` to a medium or an optical/electrical dict.

    Each copy is named after its role as well as its material ("Cu (m1)", "Cu (m2)"): Tidy3D
    warns when different media share one name.
    """
    spec = viz(role)

    def copy(m):
        return m.updated_copy(viz_spec=spec, name=f"{m.name} ({role})" if m.name else role)

    if isinstance(medium, dict):
        return {k: copy(v) for k, v in medium.items()}
    return copy(medium)
