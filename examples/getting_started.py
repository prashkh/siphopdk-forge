# %% [markdown]
# # Getting started with SiPhoPDK
#
# SiPhoPDK is a demonstration PDK for PhotonForge: a synthetic 220 nm silicon-on-insulator
# process for the O-band (1260 nm to 1360 nm), with silicon nitride, germanium photodiodes, a
# heater and two metal layers. It is not a foundry PDK. Do not use it for a tape-out.
#
# In this notebook, you will learn how to:
#
# - install SiPhoPDK,
# - check that PhotonForge finds the PDK,
# - show the technology and the cells,
# - make a Mach-Zehnder interferometer (MZI) from SiPhoPDK cells and simulate it,
# - run a design rule check (DRC) and a layout-versus-schematic check (LVS),
# - open the MZI in the local GUI.
#
# This notebook runs no FDTD simulation and no mode solver. It needs no Tidy3D account and uses
# no FlexCredits. DRC and LVS need the KLayout application, version 0.29.12 or newer
# (https://www.klayout.de/build.html).
#
# Run each cell in order: click the cell, then press **Shift + Enter**.

# %% [markdown]
# ## 1. Install SiPhoPDK
#
# This command installs SiPhoPDK from GitHub. It also installs PhotonForge 2.0, if you do not
# have it. If you imported PhotonForge in this notebook before the install, restart the kernel
# after the install (menu **Kernel**, then **Restart Kernel**).

# %%
# %pip install -q "siphopdk-forge @ git+https://github.com/prashkh/siphopdk-forge"

# If you have the wheel file, install that file instead:
# %pip install -q siphopdk_forge-1.0.0-py3-none-any.whl

# %% [markdown]
# ## 2. Check the installation
#
# A PDK tells PhotonForge that it is a PDK with an entry point in the group `photonforge.pdk`.
# The component library of the PhotonForge GUI reads this list.

# %%
from importlib.metadata import entry_points

import numpy as np
import photonforge as pf

import siphopdk

print("PhotonForge", pf.__version__)
print("SiPhoPDK", siphopdk.__version__)

pdks = {entry.name: entry.value for entry in entry_points(group="photonforge.pdk")}
print("PDKs that PhotonForge finds:", pdks)
assert "SiPhoPDK" in pdks, "PhotonForge does not find SiPhoPDK"

# %% [markdown]
# ## 3. Set the technology
#
# The technology holds the layers, the materials, the layer stack and the waveguide types (the
# port specifications). Set it as the default, so that every new component uses it.

# %%
technology = siphopdk.siphopdk()
pf.config.default_technology = technology

for name, spec in technology.ports.items():
    print(f"{name:14s} {spec.description}")

# %% [markdown]
# The cross section shows the stack: from left to right, a germanium photodiode, a ring
# modulator, a silicon waveguide, a silicon nitride waveguide, a grating coupler, an edge
# coupler and the metal stack.

# %%
ax = siphopdk.plot_cross_section()

# %% [markdown]
# ## 4. Look at the cells
#
# Every cell is a parametric component. Each one has a layout and a model that runs on your
# computer.

# %%
print(len(siphopdk.component_names), "cells:")
print(", ".join(siphopdk.component_names))

# %% [markdown]
# A cell with its default parameters, and the same cell with a smaller radius:

# %%
ring = siphopdk.component.micro_ring_modulator()
ring

# %%
small_ring = siphopdk.component.micro_ring_resonator(ring_radius=6.0)
small_ring

# %% [markdown]
# ## 5. Make an MZI
#
# Two Y-branches: the lower one splits the light, the upper one combines it. Two routes connect
# them. The outer route uses a larger bend radius, so it is longer, and the MZI is unbalanced.
#
# The routes use an analytic waveguide model with the effective index (`n_eff`) and the group
# index (`n_group`) of the SiPhoPDK strip waveguide at 1310 nm. `propagation_loss` is in dB/µm:
# 2e-4 dB/µm is 2 dB/cm.

# %%
waveguide_model = pf.AnalyticWaveguideModel(n_eff=2.59, n_group=4.29, propagation_loss=2e-4)
y_branch = siphopdk.component.y_branch()

mzi = pf.component_from_netlist(
    {
        "name": "MZI",
        "instances": {"lower": y_branch, "upper": {"component": y_branch, "origin": (0, 40)}},
        "routes": [
            (("lower", "P1"), ("upper", "P2"), {"radius": 10, "model": waveguide_model}),
            (("lower", "P2"), ("upper", "P1"), {"radius": 20, "model": waveguide_model}),
        ],
        "ports": [("lower", "P0"), ("upper", "P0")],
        "models": [pf.CircuitModel()],
    }
)
mzi

# %% [markdown]
# ## 6. Simulate the MZI
#
# The circuit model calculates the S matrix of the MZI from the models of its parts.
# `pf.C_0` is the speed of light in µm/s, so `pf.C_0 / wavelengths` gives the frequencies.

# %%
wavelengths = np.linspace(1.26, 1.36, 1001)  # µm
s_matrix = mzi.s_matrix(pf.C_0 / wavelengths)
fig, ax = pf.plot_s_matrix(s_matrix, input_ports=["P0"], output_ports=["P1"], y="dB")

# %% [markdown]
# The transmission goes up and down with the wavelength, because the two arms have different
# lengths.

# %% [markdown]
# ## 7. Check the design rules (DRC)
#
# The SiPhoPDK DRC deck has rules for width, space, enclosure, placement and area. Rule FP.1
# requires all geometry inside the floorplan layer `FLOORPLAN`, so put the MZI on a die with a
# floorplan first.

# %%
def on_a_die(component, margin=10.0):
    """Place a component on a die whose FLOORPLAN rectangle is `margin` µm larger."""
    die = pf.Component(f"{component.name}_die")
    die.add_reference(component)
    (x0, y0), (x1, y1) = component.bounds()
    die.add("FLOORPLAN", pf.Rectangle((x0 - margin, y0 - margin), (x1 + margin, y1 + margin)))
    return die


try:
    results = siphopdk.drc_check(on_a_die(mzi))
    print("MZI DRC clean:", results.is_clean)
    klayout_found = True
except (RuntimeError, OSError) as error:
    print("DRC did not run. It needs KLayout 0.29.12 or newer. The error:", error)
    klayout_found = False

# %% [markdown]
# A directional coupler with a 120 nm gap breaks rule SI.S.1: the minimum space between two
# silicon shapes is 150 nm.

# %%
if klayout_found:
    narrow = siphopdk.component.directional_coupler(coupling_distance=0.12)
    results = siphopdk.drc_check(on_a_die(narrow))
    print("Violations:", {rule: n for rule, n in results.violation_counts.items() if n})

# %% [markdown]
# ## 8. Compare the layout with the schematic (LVS)
#
# PhotonForge writes the schematic from the connections of the MZI. The LVS runset extracts the
# nets from the layout and compares the two.

# %%
if klayout_found:
    results = siphopdk.lvs_check(mzi)
    print("MZI LVS:", results.verdict)

# %% [markdown]
# LVS finds a fault that the layout does not show. Here the layout has only one of the two
# arms, and the schematic of the full MZI is the design intent:

# %%
if klayout_found:
    one_arm = pf.component_from_netlist(
        {
            "name": "MZI",
            "instances": {"lower": y_branch, "upper": {"component": y_branch, "origin": (0, 40)}},
            "routes": [(("lower", "P1"), ("upper", "P2"), {"radius": 10, "model": waveguide_model})],
            "ports": [("lower", "P0"), ("upper", "P0")],
        }
    )
    results = siphopdk.lvs_check(one_arm, schematic=mzi)
    print("One arm only, LVS:", results.verdict)

# %% [markdown]
# ## 9. Open the MZI in the local GUI
#
# The local GUI runs on this computer. Open the address that the cell prints. In the GUI, the
# component library lists SiPhoPDK with all its cells.

# %%
from photonforge.live_viewer import LiveViewer

viewer = LiveViewer(port=8765)
viewer(mzi)

# %% [markdown]
# When you are done, stop the GUI:

# %%
viewer.stop()
