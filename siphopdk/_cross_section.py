"""Stack-up visualisation for SiPhoPDK.

Builds one component holding a slice of every device family in the PDK, laid out left to
right in the order of the platform figure (Ge photodetector, microring modulator, Si WG,
SiN WG, grating coupler, edge coupler, BEOL stack), then plots the extruded cross-section.
This doubles as the acceptance check for the extrusion recipe: if a layer is missing from
the plot, its extrusion spec is wrong.
"""

import photonforge as pf

from ._viz import IMPLANT_COLORS

# x windows for each device family
_PD = (0.0, 9.0)
_MRM = (13.0, 24.0)
_SIWG = (28.0, 33.0)
_SINWG = (35.0, 38.0)
_GC = (41.0, 49.0)
_EC = (52.0, 59.0)
_BEOL = (62.0, 72.0)

_HY = 2.0  # half-height in y; the cross-section is taken at y = 0


def _rect(x0, x1, y_half=_HY):
    return pf.Rectangle((x0, -y_half), (x1, y_half))


def cross_section_component(technology=None):
    """Component whose y=0 cross-section shows the whole SiPhoPDK stack."""
    if technology is None:
        from .technology import siphopdk

        technology = siphopdk()

    c = pf.Component("SiPhoPDK stack", technology)

    # ---- Ge photodetector: slab + Ge on top + p/n contacts to M1 --------------------
    x0, x1 = _PD
    c.add(
        "SLAB", _rect(x0, x1),
        "SI", _rect(x0, x0 + 1.2), _rect(x1 - 1.2, x1),          # raised contact shelves
        "GE", _rect(x0 + 3.6, x0 + 5.4),                          # Ge absorber
        "PW", _rect(x0, x0 + 4.0), "NW", _rect(x0 + 5.0, x1),
        "PPP", _rect(x0, x0 + 1.4), "NPP", _rect(x1 - 1.4, x1),
        "VIA1", _rect(x0 + 0.3, x0 + 0.9), _rect(x1 - 0.9, x1 - 0.3),
        "M1", _rect(x0 - 0.2, x0 + 1.4), _rect(x1 - 1.4, x1 + 0.2),
    )

    # ---- Microring modulator: rib waveguide, Z-shape junction, heater above ---------
    x0, x1 = _MRM
    mid = 0.5 * (x0 + x1)
    c.add(
        "SLAB", _rect(x0, x1),
        "SI", _rect(mid - 0.25, mid + 0.25),                      # 500 nm rib core
        "SI", _rect(x0, x0 + 1.0), _rect(x1 - 1.0, x1),           # contact shelves
        "PW", _rect(x0, mid + 0.12), "NW", _rect(mid + 0.12, x1),
        "PP", _rect(x0, mid - 0.25), "NP", _rect(mid + 0.25, x1),
        "PPP", _rect(x0, x0 + 1.2), "NPP", _rect(x1 - 1.2, x1),
        "VIA1", _rect(x0 + 0.2, x0 + 0.8), _rect(x1 - 0.8, x1 - 0.2),
        "M1", _rect(x0 - 0.2, x0 + 1.2), _rect(x1 - 1.2, x1 + 0.2),
        "HTR", _rect(mid - 2.5, mid + 2.5),                       # TiN heater over the ring
        "VIA2", _rect(mid - 2.4, mid - 1.9), _rect(mid + 1.9, mid + 2.4),
        "M1", _rect(mid - 2.6, mid - 1.7), _rect(mid + 1.7, mid + 2.6),   # heater routing
    )

    # ---- passive waveguides --------------------------------------------------------
    x0, x1 = _SIWG
    c.add("SI", _rect(x0 + 0.8, x0 + 1.2))                        # 400 nm strip
    c.add("SI", _rect(x0 + 3.3, x0 + 3.8), "SLAB", _rect(x0 + 2.3, x0 + 4.8))  # 500 nm rib
    x0, x1 = _SINWG
    c.add("SIN", _rect(x0 + 1.0, x0 + 2.0))

    # ---- grating coupler: shallow-etch teeth in the silicon -------------------------
    x0, x1 = _GC
    c.add("SI", _rect(x0, x1))
    for i in range(7):
        xs = x0 + 0.7 + i * 1.0
        c.add("GRAT", _rect(xs, xs + 0.5))

    # ---- edge coupler: inverse taper tip + deep trench facet -----------------------
    x0, x1 = _EC
    c.add(
        "SI", _rect(x0 + 2.0, x1),
        "SIN", _rect(x0 + 1.0, x1),                                # SiN overlay mode expander
        "DT", _rect(x0 - 1.0, x0 + 0.5),                           # facet trench
    )

    # ---- BEOL: full metal stack and a pad opening ----------------------------------
    x0, x1 = _BEOL
    c.add(
        "M1", _rect(x0, x1),
        "VIA3", _rect(x0 + 1.0, x0 + 2.0), _rect(x1 - 2.0, x1 - 1.0),
        "M2", _rect(x0, x1),
        "PAD", _rect(x0 + 3.0, x1 - 3.0),
    )
    return c


def plot_cross_section(technology=None, ax=None, **kwargs):
    """Plot the SiPhoPDK stack-up (y = 0 cross-section of every device family)."""
    c = cross_section_component(technology)
    ax = pf.tidy3d_plot(c, y=0, ax=ax, **kwargs)
    ax.set(title=c.technology.name + " " + c.technology.version)
    return ax


_PANELS = (
    ("Ge photodetector", _PD, (-0.35, 2.4)),
    ("Microring modulator + heater", _MRM, (-0.35, 2.4)),
    ("MRM junction (zoom)", (_MRM[0] + 4.4, _MRM[0] + 7.4), (-0.12, 0.34)),
    ("Si WG (strip and rib)", _SIWG, (-0.25, 0.65)),
    ("SiN WG", _SINWG, (-0.25, 0.95)),
    ("Grating coupler", _GC, (-0.25, 0.65)),
    ("Edge coupler", _EC, (-0.9, 1.2)),
    ("BEOL: M1 / VIA3 / M2 / pad", _BEOL, (1.2, 4.6)),
)


def _silicon_z_limits(technology):
    """(slab_top, core_top) read back from the SI and SLAB extrusion specs."""
    si_layer = tuple(technology.layers["SI"].layer)
    slab_layer = tuple(technology.layers["SLAB"].layer)
    core_top = slab_top = None
    for spec in technology.extrusion_specs:
        layers = tuple(tuple(l) for l in spec.mask_layers())
        if layers == (si_layer,):
            core_top = float(spec.limits[1])
        elif layers == (slab_layer,):
            slab_top = float(spec.limits[1])
    return slab_top, core_top


def _overlay_implants(ax, component):
    """Shade the drawn implant windows, clipped to the silicon that actually exists.

    Implants carry no extrusion (they are drawing layers that feed the CHARGE doping), so
    they are invisible in the material cross-section. This overlay puts them back, which is
    also the honest way to show that the two descriptions are separate.

    An implant only exists inside silicon, so each window is split against the geometry: the
    part over the full-etch core is drawn up to the core top, the part over the partial-etch
    slab only up to the slab top. Bounding boxes are used for the drawn extents, which is
    exact for the rectangular windows a cross-section showcase uses.
    """
    from matplotlib.patches import Rectangle as MplRect

    tech = component.technology
    slab_top, core_top = _silicon_z_limits(tech)
    si = component.structures.get(tuple(tech.layers["SI"].layer), [])
    slab = component.structures.get(tuple(tech.layers["SLAB"].layer), [])
    if not (si or slab):
        return

    for name, color in IMPLANT_COLORS.items():
        windows = component.structures.get(tuple(tech.layers[name].layer), [])
        if not windows:
            continue
        in_core = pf.boolean(windows, si, "*") if si else []
        on_slab = pf.boolean(pf.boolean(windows, slab, "*"), si, "-") if slab else []
        for shapes, z_top in ((in_core, core_top), (on_slab, slab_top)):
            if z_top is None:
                continue
            for shape in shapes:
                (x0, _), (x1, _) = shape.bounds()
                ax.add_patch(
                    MplRect((x0, 0.0), x1 - x0, z_top, facecolor=color, edgecolor=color,
                            linewidth=0.6, alpha=0.42, zorder=5)
                )


def plot_stack_panels(technology=None, figsize=(15, 8), show_implants=True):
    """Zoomed per-device panels of the stack, in the style of a platform cross-section figure."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    c = cross_section_component(technology)
    n = len(_PANELS)
    rows = (n + 2) // 2 if show_implants else (n + 1) // 2  # spare cell for the legend
    fig, axes = plt.subplots(rows, 2, figsize=figsize)
    axes = axes.ravel()
    for ax, (title, (x0, x1), zlim) in zip(axes, _PANELS):
        pf.tidy3d_plot(c, y=0, ax=ax)
        if show_implants:
            _overlay_implants(ax, c)
        pad = 0.6 if (x1 - x0) > 4 else 0.0
        ax.set(xlim=(x0 - pad, x1 + pad), ylim=zlim, title=title, xlabel="x (um)", ylabel="z (um)")
    for ax in axes[n:]:
        ax.axis("off")
    if show_implants:
        handles = [Patch(facecolor=col, alpha=0.6, label=k) for k, col in IMPLANT_COLORS.items()]
        axes[n].legend(handles=handles, loc="center", ncol=3, frameon=False,
                       title="implants (drawing layers, not extruded)")
    fig.suptitle(f"{c.technology.name} {c.technology.version} stack-up", fontsize=13)
    fig.tight_layout()
    return fig


def layer_coverage(technology=None):
    """Report which layers are structural (have an extrusion) and which are not."""
    if technology is None:
        from .technology import siphopdk

        technology = siphopdk()
    used = set()
    for spec in technology.extrusion_specs:
        for lay in spec.mask_layers():
            used.add(tuple(lay))
    rows = []
    for name, spec in technology.layers.items():
        rows.append((name, tuple(spec.layer), tuple(spec.layer) in used))
    return rows
