"""Compute the mode indices of the SiPhoPDK waveguides, and check the ring calibration.

The ring modulator is calibrated to Y. Yuan et al., Nat. Commun. 15, 918 (2024): a 12 um radius
ring in a 500 nm rib waveguide (220 nm Si, 90 nm slab) at 1310 nm has FSR = 5.7 nm. With
FSR = lambda^2 / (n_g * 2 * pi * R), the rib mode then needs a group index of

    n_g = 1.310^2 / (0.0057 * 2 * pi * 12) = 3.99

The values this script prints are the ones in ``siphopdk/_calibration.py``.

    python scripts/calibrate.py

CAUTION: the script uses PhotonForge's local mode solver, which needs the ``tidy3d-extras``
package. Without it, the mode solves run in the Flexcompute cloud and use FlexCredits.
"""

import numpy as np
import photonforge as pf

import siphopdk

pf.config.use_local_mode_solver = True

LAMBDA = 1.310
RADIUS = 12.0
FSR_TARGET = 0.0057
NG_TARGET = LAMBDA**2 / (FSR_TARGET * 2 * np.pi * RADIUS)

tech = siphopdk.siphopdk()
pf.config.default_technology = tech
freqs = pf.C_0 / np.array([1.290, 1.310, 1.330])

print(f"target n_g for an FSR of {FSR_TARGET * 1e3:.1f} nm at R = {RADIUS} um: {NG_TARGET:.3f}\n")

rows = []
for name in ("SI_RIB_O", "SI_STRIP_O", "SIN_STRIP_O"):
    solver = pf.port_modes(
        tech.ports[name], freqs, group_index=True, mesh_refinement=25, verbose=False,
        show_progress=False,
    )
    n_eff = np.real(solver.data.n_eff.values[:, 0])
    n_g = np.real(solver.data.n_group.values[:, 0])
    fsr = LAMBDA**2 / (n_g[1] * 2 * np.pi * RADIUS)
    rows.append((name, n_eff[1], n_g[1], fsr))
    print(f"{name:14s} n_eff = {n_eff[1]:.4f}   n_g = {n_g[1]:.4f}   ring FSR = {fsr * 1e3:.2f} nm")

rib = rows[0]
error = 100 * (rib[3] - FSR_TARGET) / FSR_TARGET
print(f"\nrib: n_g {rib[2]:.3f} against {NG_TARGET:.3f}, FSR {rib[3] * 1e3:.2f} nm against 5.70 nm "
      f"({error:+.1f} %)")
