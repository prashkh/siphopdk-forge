"""SiPhoPDK calibration constants.

Optical values are computed on SiPhoPDK's own stack with the Tidy3D mode solver
(``scripts/calibrate.py`` in the source repository). Device values are anchored to Yuan et
al., Nat. Commun. 15:918 (2024), the 5 x 200 Gbps two-segment Z-junction MRM chip this PDK is
built to reproduce.

Nothing here is fitted to a foundry PDK.
"""

import numpy as np
import photonforge as pf

WAVELENGTH = 1.310  # um, O-band design point
REF_FREQUENCY = pf.C_0 / WAVELENGTH

# --- waveguides, from scripts/calibrate.py on this stack --------------------------
N_EFF_STRIP, N_GROUP_STRIP = 2.5932, 4.2949  # SI_STRIP_O, 400 nm
N_EFF_RIB, N_GROUP_RIB = 2.7856, 3.9854      # SI_RIB_O, 500 nm on a 90 nm slab
N_EFF_SIN, N_GROUP_SIN = 1.6429, 2.0302      # SIN_STRIP_O, 1 um

# [SYNTH] propagation loss, dB/um. 2 dB/cm strip, 0.5 dB/cm rib-mode SiN.
LOSS_STRIP = 2.0e-4
LOSS_RIB = 2.5e-4
LOSS_SIN = 0.5e-4
# [SYNTH] excess loss per 90 degree bend and per crossing, dB
LOSS_BEND_90 = 0.01
LOSS_CROSSING = 0.13
CROSSTALK_CROSSING = 10 ** (-40 / 20)

# --- ring modulator, anchored to Yuan 2024 ----------------------------------------
R_RING = 12.0
L_RING = 2 * np.pi * R_RING
FSR_NM = 5.7
N_GROUP_RING = WAVELENGTH**2 / (FSR_NM * 1e-3 * L_RING)  # 3.99, matches the rib mode
ALPHA0_DBUM = 0.0129   # 129 dB/cm inside the doped ring (paper SI)
KAPPA1 = 0.387         # bus-ring power coupling, from the paper's Q and ER
DLAMBDA_DV_PM = 16.3   # MSB segment, paper Fig. 2h
DN_DV = (N_GROUP_RING / WAVELENGTH) * DLAMBDA_DV_PM * 1e-6
F_3DB_RC = 64e9        # MSB equivalent circuit, paper Fig. 3e
V_PI_L = 0.53e4        # V.um, 0.53 V.cm for the MSB Z-shape junction

# --- WDM channel plan ---------------------------------------------------------------
# A ring's resonance moves with its circumference, but NOT as lambda/L: the resonance
# condition n_eff(lambda) L = m lambda differentiates to
#     d(lambda)/dL = (lambda / L) * (n_eff / n_g)
# so the radius step for a wanted channel spacing carries the dispersion factor n_g/n_eff
# (1.43 on this stack). Using the naive lambda/L gives 0.79 nm where 1.14 nm was intended.
CHANNEL_SPACING_HZ = 200e9
CHANNEL_SPACING_NM = WAVELENGTH**2 * CHANNEL_SPACING_HZ / pf.C_0 * 1e3  # 1.145 nm
NUM_CHANNELS = 5
R_CH0 = 11.9942          # puts a resonance at exactly 1310.0 nm on this stack
DR_PER_CHANNEL = R_CH0 * (CHANNEL_SPACING_NM * 1e-3 / WAVELENGTH) * (N_GROUP_RIB / N_EFF_RIB)
CHANNEL_RADII = [round(R_CH0 + i * DR_PER_CHANNEL, 4) for i in range(NUM_CHANNELS)]
# The grid those radii are aiming at, in um. Derived rather than written down, so the layout
# and anything that simulates it cannot disagree about where a channel is supposed to be: the
# radii come from the spacing and so does this. Measured on the built cells at 1310.002,
# 1311.147, 1312.290, 1313.433 and 1314.575 nm, ER 15.8-16.0 dB.
CHANNEL_WAVELENGTHS = [WAVELENGTH + i * CHANNEL_SPACING_NM * 1e-3 for i in range(NUM_CHANNELS)]

# --- thermal ----------------------------------------------------------------------
DN_DT = 1.86e-4        # /K
R_TH = 1.20            # K/mW, from the paper's 73 pm/mW heater

# --- germanium photodiode [PAPER: platform paper Fig. 8] --------------------------
PD_RESPONSIVITY = 0.95   # A/W, paper reports close to 1 A/W
PD_DARK_CURRENT = 15e-9  # A, paper reports < 20 nA
PD_BANDWIDTH = 110e9     # Hz, paper reports 110 GHz at 200 uA photocurrent

# --- monitor tap [SYNTH] ----------------------------------------------------------
TAP_RATIO = 0.01         # 1 % to the monitor arm

# --- grating coupler [SYNTH], shaped like the paper's SPGC ------------------------
GC_PEAK_WAVELENGTH = 1.310
GC_PEAK_LOSS_DB = 1.3   # paper Fig. 6c: peak loss < 1.3 dB
GC_BW_1DB_NM = 20.0     # paper Fig. 6c: 1 dB bandwidth about 20 nm

# --- Tier 2 ------------------------------------------------------------------------
# [SYNTH] Si-to-SiN interlayer transition, dB. Adiabatic over 40 um at this index contrast;
# published couplers of this kind report 0.1-0.5 dB, so this sits in the middle.
LOSS_SI_TO_SIN = 0.25
# [SYNTH] SiN inverse-taper edge coupler to a lensed fibre, dB, and its wavelength slope.
# Broadband is the whole reason to prefer one over a grating: 0.02 dB/nm across the O-band is
# under half a dB edge to edge, where the grating loses 1 dB in 20 nm.
EC_PEAK_LOSS_DB = 1.5
EC_SLOPE_DB_PER_NM = 0.02

# [PAPER] Loaded Q of the ring modulator. Yuan 2024 measures 3700, and it is not decoration:
# a cavity holds light for Q/(2*pi*f) seconds, so it is also a low-pass filter on anything
# modulated inside it, with a 3 dB point at f/Q - about 60 GHz here. At 200 Gb/s that is the
# narrowest bandwidth in the link, which is the trade a ring modulator makes for its size.
Q_LOADED = 3700.0
