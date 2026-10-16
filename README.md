# SiPhoPDK for PhotonForge

SiPhoPDK is a synthetic 220 nm silicon-on-insulator process for the O-band, with a silicon
nitride layer, germanium photodetectors, a TiN heater and two metal layers. It is a
demonstration PDK for [PhotonForge](https://docs.flexcompute.com/projects/photonforge/): it
holds the technology, 27 parametric cells with offline models, a KLayout DRC deck and a KLayout
LVS runset.

SiPhoPDK is not a foundry PDK. Every number in it is either published physics, taken from an
open source, or invented for teaching. Do not use it for a tape-out.


## Installation

Installation via `pip`:

    pip install siphopdk-forge

The package needs PhotonForge 2.0 or newer. DRC and LVS also need the
[KLayout](https://www.klayout.de/build.html) application, version 0.29.12 or newer.


## Usage

Set the SiPhoPDK technology as the default:

    import photonforge as pf
    import siphopdk

    tech = siphopdk.siphopdk()
    pf.config.default_technology = tech

The `siphopdk` function creates a parametric technology. Its parameters set the layer
thicknesses, the etch sidewall angle and the mask bias, for example
`siphopdk.siphopdk(si_thickness=0.218)`.

The cells are in the `component` submodule:

    print(siphopdk.component_names)

    ring = siphopdk.component.micro_ring_modulator(ring_radius=12.0)

Each function has its own documentation:

    help(siphopdk.siphopdk)

    help(siphopdk.component.micro_ring_modulator)

To show the extruded cross section of the stack:

    siphopdk.plot_cross_section()

When the package is installed, the component library of the PhotonForge GUI lists SiPhoPDK
with all its cells.

`examples/getting_started.ipynb` installs SiPhoPDK, shows the technology and the cells,
makes and simulates an MZI, runs DRC and LVS on it and opens it in the GUI.


## Port specifications

| Name | Waveguide |
|---|---|
| `SI_STRIP_O` | Si strip, 400 nm, TE |
| `SI_STRIP_O_TM` | Si strip, 400 nm, TM |
| `SI_RIB_O` | Si rib, 500 nm on a 90 nm slab, TE |
| `SI_WIDE_O` | Si strip, 1000 nm, fundamental TE mode |
| `SI_MM_O` | Si multimode strip, 2000 nm |
| `SIN_STRIP_O` | SiN strip, 1000 nm, TE |


## Cells

Every cell is a parametric component with a model that runs on your computer: an analytic
model, a ring model, or a data model built from synthetic measurement data. No cell needs a
simulation before you can use it in a circuit.

| Cell | Description | Ports | Terminals |
|---|---|---|---|
| `straight` | Straight waveguide | P0, P1 | |
| `bend` | Euler-assisted circular bend | P0, P1 | |
| `taper` | Linear width taper between two port specifications | P0, P1 | |
| `strip_to_rib` | Strip-to-rib transition | P0, P1 | |
| `crossing` | Waveguide crossing | P0 to P3 | |
| `y_branch` | 1x2 splitter | P0 to P2 | |
| `mmi_2x2` | 2x2 multimode interference coupler, 50:50 | P0 to P3 | |
| `directional_coupler` | Symmetric directional coupler | P0 to P3 | |
| `tap` | Monitor tap, 1 % to the drop arm | P0 to P3 | |
| `termination` | Radiative termination | P0 | |
| `grating_coupler` | Focusing grating coupler | P0, FIBER | |
| `sin_edge_coupler` | Inverse-taper edge coupler on SiN | P0, FACET | |
| `si_to_sin` | Silicon to SiN interlayer transition | P0, P1 | |
| `spiral_delay` | Double-spiral delay line of exact length | P0, P1 | |
| `contra_dc` | Contra-directional coupler (Bragg add-drop) | P0 to P3 | |
| `grating_antenna` | Weak grating emitter for a phased array | P0 | |
| `thermal_phase_shifter` | Strip waveguide with a TiN heater | P0, P1 | HTRP, HTRN |
| `pn_phase_shifter` | Rib waveguide with a depletion-mode p-n junction | P0, P1 | ANODE, CATHODE |
| `pin_attenuator` | Carrier-injection variable optical attenuator | P0, P1 | ANODE, CATHODE |
| `micro_ring_resonator` | Ring resonator with a heater | P0, P1 | HTRP, HTRN |
| `micro_ring_modulator` | Depletion-mode ring modulator with a heater | P0, P1 | ANODE, CATHODE, HTRP, HTRN |
| `add_drop_ring` | Add-drop ring with an optional heater | P0 to P3 | HTRP, HTRN |
| `add_drop_ring_corner` | Add-drop ring between a row bus and a column bus | P0 to P3 | HTRP, HTRN |
| `photodiode` | Germanium photodiode, lateral p-i-n | P0 | ANODE, CATHODE |
| `tw_electrode` | Travelling-wave GSG electrode on M2 | | SIN, SOUT, 4 grounds |
| `bond_pad` | Bond pad on M2, with an M1 landing | | M1, T |
| `via_stack` | M1 to M2 via stack | | M1, M2 |


## Design rule check

The DRC deck has 43 rules for width, space, enclosure, placement and area. Each rule in the deck
says where its value comes from: the published literature, an open PDK, the physics of the
SiPhoPDK stack, or a value chosen for teaching. Run it on a component or on a layout file:

    results = siphopdk.drc_check(chip)
    print(results.is_clean)
    print(results)

`siphopdk.DRC_RULES` is the path of the deck, `rules/siphopdk.drc`. Use it with
`photonforge.drc_check`, in the DRC panel of the PhotonForge GUI, or with KLayout directly:

    klayout -b -r siphopdk.drc -rd gdsfile=chip.gds -rd resultsfile=chip.lyrdb

Geometry inside the `EXCL_DRC` layer (95/0) is waived.


## Layout versus schematic

The LVS runset extracts the optical and the electrical nets from the layout and compares them
with a schematic. PhotonForge writes the schematic from the connections of the component, so a
circuit of SiPhoPDK cells and PhotonForge routes is checked with one call:

    results = siphopdk.lvs_check(chip)
    print(results.verdict)

To check a layout against a different design intent, give that intent as the schematic:

    results = siphopdk.lvs_check(layout, schematic=intended_circuit)

`siphopdk.LVS_RULES` is the path of the runset, `rules/siphopdk.lvs`. It reads these layers,
which every SiPhoPDK cell draws:

| Layer | Use |
|---|---|
| `DEVREC` (100/0) | Outline of each cell, which ends at its port planes |
| `PINREC` (101/0) | A 200 nm box on each optical port |
| `PINREC_E` (102/0) | A 200 nm box on each electrical terminal |
| text on 0/0 | The name of each port and terminal, on its box |

A waveguide outside every outline joins the pins it touches. Only the contact implants, the vias
and the two metals conduct: a heater and a p-n junction show as two pins, so a trace that
bridges them is reported as a short.


## Warnings

The 3D structures that PhotonForge makes from this technology approximate the structures a
process would make. Doping profiles are solids with sharp boundaries and uniform concentration.
The cell models use synthetic data where the source gives no measured value. SiPhoPDK is for
teaching and for trying out a design flow, not for predicting the performance of a real chip.


## References

- Y. Yuan et al., "A 5 x 200 Gbps microring modulator silicon chip empowered by two-segment
  Z-shape junctions", Nat. Commun. 15, 918 (2024),
  [doi:10.1038/s41467-024-45301-3](https://doi.org/10.1038/s41467-024-45301-3). The ring
  modulator, its junction and its model are calibrated to this paper.
- S. K. Yeh et al., "Silicon Photonics Platform for Next Generation Data Communication
  Technologies". The process stack follows the generic platform in this paper.
- [SiEPIC EBeam PDK](https://github.com/SiEPIC/SiEPIC_EBeam_PDK): the open rules marked
  `[OPEN]` in the DRC deck.


## Changelog

### 1.0.0 - 2026-10-15

- First public release.
- Technology, 6 port specifications and 27 parametric cells with offline models.
- KLayout DRC deck (43 rules) and LVS runset, with `drc_check` and `lvs_check`.
- Design rules on the layer specifications (`min_width`, `min_spacing`), conductive layers and
  layer connections, for the PhotonForge routers.


## License

SiPhoPDK is free software under the GNU General Public License, version 3 or later. See
[LICENSE](LICENSE) for the license text and [NOTICE](NOTICE) for the copyright notice.

Under section 7(b) of the license, a work that contains SiPhoPDK, or is based on it, must keep
the attribution "SiPhoPDK by Flexcompute Inc." in its copyright notices and in the legal notices
it displays.


## Citation

If you use SiPhoPDK in published work, please cite it as given in [CITATION.cff](CITATION.cff).
