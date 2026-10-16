"""Shared test helpers: is KLayout usable, and the DRC violations of a cell, by rule."""

import photonforge as pf
import pytest

import siphopdk

MIN_KLAYOUT = (0, 29, 12)


def _klayout_ok():
    try:
        from photonforge.extension import _find_klayout, _klayout_version

        version = _klayout_version(_find_klayout())
        numbers = tuple(int(x) for x in version.split()[-1].split(".")[:3])
        return numbers >= MIN_KLAYOUT
    except Exception:
        return False


KLAYOUT_OK = _klayout_ok()

needs_klayout = pytest.mark.skipif(
    not KLAYOUT_OK, reason="needs the KLayout application, 0.29.12 or newer"
)


def drc_counts(*items, name="drc_check"):
    """Place the cells or references in a floorplan and return {rule: count} for every rule hit.

    A plain cell is placed at the origin. The FLOORPLAN rectangle reaches 5 um past the
    content, so FP.1 holds.
    """
    technology = siphopdk.siphopdk()
    wrapper = pf.Component(name, technology)
    for item in items:
        wrapper.add(item if isinstance(item, pf.Reference) else pf.Reference(item))
    (x0, y0), (x1, y1) = wrapper.bounds()
    wrapper.add("FLOORPLAN", pf.Rectangle((x0 - 5, y0 - 5), (x1 + 5, y1 + 5)))
    results = siphopdk.drc_check(wrapper, show_progress=False)
    return {rule: n for rule, n in results.violation_counts.items() if n}
