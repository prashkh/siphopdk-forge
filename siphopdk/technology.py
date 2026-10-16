"""SiPhoPDK technology: a synthetic 220 nm SOI + SiN + Ge O-band datacom platform.

The stack follows the generic platform described in the public literature (S. K. Yeh et
al., "Silicon Photonics Platform for Next Generation Data Communication Technologies",
Fig. 1b) and the 220 nm SOI cross-section published in Yuan et al., Nat. Commun. 15:918
(2024), which is the chip this PDK is calibrated to reproduce. Nothing here is derived
from a foundry PDK.

Calibration target (Yuan 2024): a 12 um radius ring in the 500 nm rib waveguide has
FSR = 5.7 nm at 1310 nm, which requires n_g = 3.99. ``scripts/calibrate.py`` in the source
repository computes the mode indices of the port specifications with the Tidy3D mode solver.
"""

import photonforge as pf
import photonforge.typing as pft
import tidy3d as td

from ._layers import _layers
from ._viz import with_viz

_L = {k: v.layer for k, v in _layers.items()}  # name -> (layer, datatype)

_SIO2_BASE = {
    "optical": td.material_library["SiO2"]["Palik_NoLoss"],
    "electrical": td.Medium(permittivity=3.9, name="SiO2"),
}
_SI_BASE = {
    "optical": td.material_library["cSi"]["Li1993_293K"],
    "electrical": td.Medium(permittivity=11.7, name="Si"),
}
_SIN_BASE = {
    "optical": td.material_library["Si3N4"]["Luke2015Sellmeier"],
    "electrical": td.Medium(permittivity=7.0, name="Si3N4"),
}
_GE_BASE = {
    # Nunley covers 0.19 to 2.48 um and absorbs at 1310 nm (k = 0.086); Icenogle1976 starts at
    # 2.5 um and would have no absorption in the O-band.
    "optical": td.material_library["Ge"]["Nunley"],
    "electrical": td.Medium(permittivity=16.0, name="Ge"),
}
# TiN has no entry in the Tidy3D material library; W is the closest optical stand-in and
# the electrical model is what matters for the heater and RF work.
_TIN_BASE = {
    "optical": td.material_library["W"]["Werner2009"],
    "electrical": td.LossyMetalMedium(
        conductivity=4.5,  # S/um, TiN ~ 22 uOhm.cm
        frequency_range=[0.1e9, 200e9],
        fit_param=td.SurfaceImpedanceFitterParam(max_num_poles=16),
    ),
}
_CU_BASE = {
    "optical": td.material_library["Cu"]["RakicLorentzDrude1998"],
    "electrical": td.LossyMetalMedium(
        conductivity=58.0,  # S/um
        frequency_range=[0.1e9, 200e9],
        fit_param=td.SurfaceImpedanceFitterParam(max_num_poles=16),
    ),
}
_AIR_BASE = td.Medium(permittivity=1.0)

# One medium instance per visual role: Tidy3D colours structures by medium, so M1 and M2
# (both copper) must be distinct objects to plot in different colours.
_SUBSTRATE = with_viz(_SI_BASE, "substrate")
_OXIDE = with_viz(_SIO2_BASE, "oxide")
_SI_CORE = with_viz(_SI_BASE, "si_core")
_SI_SLAB = with_viz(_SI_BASE, "si_slab")
_SIN = with_viz(_SIN_BASE, "sin")
_GE = with_viz(_GE_BASE, "ge")
_HEATER = with_viz(_TIN_BASE, "heater")
_CONTACT = with_viz(_CU_BASE, "contact")
_M1 = with_viz(_CU_BASE, "m1")
_VIA = with_viz(_CU_BASE, "via")
_M2 = with_viz(_CU_BASE, "m2")
_OPENING = with_viz(_AIR_BASE, "opening")
_TRENCH = with_viz(_AIR_BASE, "trench")


@pf.parametric_technology
def siphopdk(
    *,
    si_thickness: pft.PositiveDimension = 0.220,
    slab_thickness: pft.PositiveDimension = 0.090,
    grating_etch_depth: pft.PositiveDimension = 0.070,
    sin_thickness: pft.PositiveDimension = 0.300,
    sin_gap: pft.PositiveDimension = 0.100,
    ge_thickness: pft.PositiveDimension = 0.250,
    heater_thickness: pft.PositiveDimension = 0.120,
    heater_gap: pft.PositiveDimension = 0.400,
    m1_thickness: pft.PositiveDimension = 0.500,
    m1_gap: pft.PositiveDimension = 0.460,
    m2_thickness: pft.PositiveDimension = 1.000,
    m2_gap: pft.PositiveDimension = 0.800,
    passivation_thickness: pft.PositiveDimension = 0.500,
    box_thickness: pft.PositiveDimension = 2.000,
    sidewall_angle: pft.Angle = 5.0,
    si_mask_dilation: pft.Coordinate = 0.0,
    slab_mask_dilation: pft.Coordinate = 0.0,
    sin_mask_dilation: pft.Coordinate = 0.0,
    rib_width: pft.PositiveDimension = 0.500,
    rib_slab_width: pft.PositiveDimension = 3.000,
    strip_width: pft.PositiveDimension = 0.400,
    sin_width: pft.PositiveDimension = 1.000,
) -> pf.Technology:
    """SiPhoPDK: synthetic 220 nm SOI + SiN O-band platform.

    Args:
        si_thickness: Full silicon core thickness.
        slab_thickness: Silicon thickness left in the partial-etch (slab) region.
        grating_etch_depth: Shallow etch depth used by grating couplers.
        sin_thickness: Silicon nitride thickness.
        sin_gap: Oxide gap between the top of the silicon and the bottom of the nitride.
        ge_thickness: Germanium epi thickness, grown on the slab (photodetectors).
        heater_thickness: TiN heater thickness.
        heater_gap: Oxide gap between the top of the nitride and the heater.
        m1_thickness: Metal 1 thickness.
        m1_gap: Oxide gap between the top of the heater and Metal 1.
        m2_thickness: Metal 2 thickness (RF and power).
        m2_gap: Oxide gap between the top of Metal 1 and Metal 2.
        passivation_thickness: Oxide above Metal 2.
        box_thickness: Buried oxide thickness.
        sidewall_angle: Etch sidewall angle for silicon and nitride, in degrees.
        si_mask_dilation: Bias applied to the full-etch silicon mask (CD control).
        slab_mask_dilation: Bias applied to the slab mask.
        sin_mask_dilation: Bias applied to the nitride mask.
        rib_width: Width of the standard rib waveguide (modulators).
        rib_slab_width: Total slab width of the standard rib waveguide.
        strip_width: Width of the standard single-mode strip waveguide.
        sin_width: Width of the standard nitride waveguide.

    Returns:
        Technology: SiPhoPDK technology definition.
    """
    layers = {k: v.copy() for k, v in _layers.items()}

    z_si_top = si_thickness
    z_ge_bot = slab_thickness  # Ge is grown in a window on the slab
    z_ge_top = z_ge_bot + ge_thickness
    z_sin = z_si_top + sin_gap
    z_sin_top = z_sin + sin_thickness
    z_htr = z_sin_top + heater_gap
    z_htr_top = z_htr + heater_thickness
    z_m1 = z_htr_top + m1_gap
    z_m1_top = z_m1 + m1_thickness
    z_m2 = z_m1_top + m2_gap
    z_m2_top = z_m2 + m2_thickness
    z_top = z_m2_top + passivation_thickness

    extrusion_specs = [
        # substrate and cladding
        pf.ExtrusionSpec(pf.MaskSpec(), _SUBSTRATE, (-pf.Z_INF, -box_thickness)),
        pf.ExtrusionSpec(pf.MaskSpec(), _OXIDE, (-box_thickness, z_top)),
        # silicon: full etch, then the partial-etch slab
        pf.ExtrusionSpec(
            pf.MaskSpec(_L["SI"], dilation=si_mask_dilation), _SI_CORE, (0, si_thickness), sidewall_angle
        ),
        pf.ExtrusionSpec(
            pf.MaskSpec(_L["SLAB"], dilation=slab_mask_dilation),
            _SI_SLAB,
            (0, slab_thickness),
            sidewall_angle,
        ),
        # grating shallow etch removes silicon from the top down
        pf.ExtrusionSpec(
            pf.MaskSpec(_L["GRAT"]), _OXIDE, (si_thickness - grating_etch_depth, si_thickness)
        ),
        # germanium on top of the silicon
        pf.ExtrusionSpec(pf.MaskSpec(_L["GE"]), _GE, (z_ge_bot, z_ge_top)),
        # nitride
        pf.ExtrusionSpec(
            pf.MaskSpec(_L["SIN"], dilation=sin_mask_dilation),
            _SIN,
            (z_sin, z_sin_top),
            sidewall_angle,
        ),
        # heater, vias, metals
        pf.ExtrusionSpec(pf.MaskSpec(_L["HTR"]), _HEATER, (z_htr, z_htr_top)),
        pf.ExtrusionSpec(pf.MaskSpec(_L["VIA1"]), _CONTACT, (slab_thickness, z_m1)),
        pf.ExtrusionSpec(pf.MaskSpec(_L["VIA2"]), _CONTACT, (z_htr_top, z_m1)),
        pf.ExtrusionSpec(pf.MaskSpec(_L["M1"]), _M1, (z_m1, z_m1_top)),
        pf.ExtrusionSpec(pf.MaskSpec(_L["VIA3"]), _VIA, (z_m1_top, z_m2)),
        pf.ExtrusionSpec(pf.MaskSpec(_L["M2"]), _M2, (z_m2, z_m2_top)),
        # pad opening and deep trench
        pf.ExtrusionSpec(pf.MaskSpec(_L["PAD"]), _OPENING, (z_m2_top, pf.Z_INF)),
        pf.ExtrusionSpec(pf.MaskSpec(_L["DT"]), _TRENCH, (-pf.Z_INF, pf.Z_INF)),
    ]

    _opt_lo, _opt_hi = -0.9, si_thickness + 0.9
    ports = {
        "SI_STRIP_O": pf.PortSpec(
            description=f"Si strip TE, O-band, w={strip_width * 1e3:.0f} nm",
            width=2.4,
            limits=(_opt_lo, _opt_hi),
            num_modes=1,
            polarization=None,
            target_neff=3.5,
            path_profiles=((strip_width, 0.0, _L["SI"]),),
            default_radius=10.0,   # the routing default: a 400 nm strip is tightly confined, and the
            #                 PDK's own bend cell is characterised here
        ),
        "SI_STRIP_O_TM": pf.PortSpec(
            description=f"Si strip TM, O-band, w={strip_width * 1e3:.0f} nm",
            width=2.4,
            limits=(_opt_lo, _opt_hi),
            num_modes=1,
            added_solver_modes=1,
            polarization="TM",
            target_neff=3.5,
            path_profiles=((strip_width, 0.0, _L["SI"]),),
            default_radius=15.0,   # TM is less confined than TE and leaks sooner in a bend
        ),
        "SI_RIB_O": pf.PortSpec(
            description=f"Si rib TE, O-band, w={rib_width * 1e3:.0f} nm on {slab_thickness * 1e3:.0f} nm slab",
            width=rib_slab_width + 2.0,
            limits=(_opt_lo, _opt_hi),
            num_modes=1,
            polarization=None,
            target_neff=3.5,
            path_profiles=((rib_width, 0.0, _L["SI"]), (rib_slab_width, 0.0, _L["SLAB"])),
            default_radius=15.0,   # the slab carries light sideways, so a rib bends worse than a strip
        ),
        "SI_WIDE_O": pf.PortSpec(
            description="Si wide strip TE, O-band, w=1000 nm, fundamental mode only",
            width=3.5,
            limits=(_opt_lo, _opt_hi),
            num_modes=1,
            polarization=None,
            target_neff=3.5,
            path_profiles=((1.0, 0.0, _L["SI"]),),
            default_radius=10.0,   # wider than the strip and no harder to bend
        ),
        "SI_MM_O": pf.PortSpec(
            description="Si multimode strip TE, O-band, w=2000 nm",
            width=6.0,
            limits=(-1.5, si_thickness + 1.5),
            num_modes=8,
            polarization=None,
            target_neff=3.5,
            path_profiles=((2.0, 0.0, _L["SI"]),),
            default_radius=50.0,   # a multimode guide bent tightly converts power between its modes
        ),
        "SIN_STRIP_O": pf.PortSpec(
            description=f"SiN strip TE, O-band, w={sin_width * 1e3:.0f} nm",
            width=4.0,
            limits=(z_sin - 1.2, z_sin_top + 1.2),
            num_modes=1,
            polarization=None,
            target_neff=2.0,
            path_profiles=((sin_width, 0.0, _L["SIN"]),),
            default_radius=60.0,   # nitride's index contrast is a third of silicon's, so the
            #                 mode is loosely held and a tight bend simply radiates
        ),
    }

    # Layer pairs whose geometry joins electrically. The contact implants reach M1 through
    # VIA1, the heater reaches M1 through VIA2 (the heater film itself is a resistor and is not
    # in this table), and M1 reaches M2 through VIA3.
    connected_layers = [
        (_L["NPP"], _L["VIA1"]),
        (_L["PPP"], _L["VIA1"]),
        (_L["VIA1"], _L["M1"]),
        (_L["VIA2"], _L["M1"]),
        (_L["M1"], _L["VIA3"]),
        (_L["VIA3"], _L["M2"]),
    ]

    return pf.Technology(
        "SiPhoPDK",
        "1.0.0",
        layers,
        extrusion_specs,
        ports,
        background_medium=_OXIDE,
        connected_layers=connected_layers,
    )
