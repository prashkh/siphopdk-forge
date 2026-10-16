"""The DRC deck in `siphopdk/rules/siphopdk.drc`, run through `siphopdk.drc_check`.

Four checks: every rule fires on its gallery site, a compliant cell is clean, EXCL_DRC waives
what it covers and nothing else, and every width, space and enclosure rule fires exactly at its
limit (1 nm under it is flagged, at it is not). Then every PDK cell, with its defaults, must be
clean.
"""

import re

import photonforge as pf
import pytest
from drc_gallery import BUILDERS, build_clean, build_gallery, build_waiver
from helpers import drc_counts, needs_klayout

import siphopdk

pytestmark = needs_klayout

DECK = siphopdk.DRC_RULES.read_text()
RULE_IDS = re.findall(r'\.output\("([^"]+)"', DECK)


@pytest.fixture(scope="module")
def technology():
    tech = siphopdk.siphopdk()
    pf.config.default_technology = tech
    return tech


def check(component):
    return siphopdk.drc_check(component, show_progress=False)


def test_deck_has_43_rules_and_a_gallery_site_for_each():
    assert len(RULE_IDS) == 43
    assert len(set(RULE_IDS)) == len(RULE_IDS), "a rule id is used twice"
    assert set(RULE_IDS) == set(BUILDERS), sorted(set(RULE_IDS) ^ set(BUILDERS))


def test_every_rule_fires_on_the_gallery(technology):
    counts = check(build_gallery(technology)).violation_counts
    silent = [rule for rule in RULE_IDS if counts.get(rule, 0) == 0]
    assert not silent, f"rules that did not fire: {silent}"


def test_the_compliant_cell_is_clean(technology):
    results = check(build_clean(technology))
    assert results.is_clean, {r: n for r, n in results.violation_counts.items() if n}


def test_excl_drc_waives_only_what_it_covers(technology):
    counts = check(build_waiver(technology)).violation_counts
    assert counts.get("SI.W.1") == 1, "one violation is waived and one is reported"


# ---------------------------------------------------------------------------------------
# Thresholds: each width, space and enclosure rule, drawn 1 nm under its limit and at it
# ---------------------------------------------------------------------------------------
SITE = 300.0  # um between test sites: beyond every rule distance and every test shape


def _layers():
    return {m[0]: (int(m[1]), int(m[2])) for m in re.findall(r"^(\w+)\s*=\s*input\((\d+),\s*(\d+)\)", DECK, re.M)}


def _cases():
    """(rule, kind, layers, limit) for every width, space and enclosure rule of the deck."""
    layers = _layers()
    cases = []
    for line in DECK.splitlines():
        rule = re.search(r'\.output\("([^"]+)"', line)
        if not rule or rule[1] == "SLAB.E.1":  # SLAB.E.1 measures a derived wing: gallery only
            continue
        width = re.search(r"waive\(\(?(\w+)[^.]*\.(width|space)\(([\d.]+)", line)
        enclosure = re.search(r"waive\((\w+)\.enclosed\((\w+),\s*([\d.]+)", line)
        if enclosure:
            cases.append((rule[1], "enclosure", (layers[enclosure[1]], layers[enclosure[2]]), float(enclosure[3])))
        elif width:
            cases.append((rule[1], width[2], (layers[width[1]],), float(width[3])))
    return cases


def _site(c, kind, layers, d, ox, oy, limit):
    """Draw one test site at (ox, oy) with distance d; return its bounding box."""
    length = max(5.0, 3 * limit)
    if kind == "width":
        c.add(layers[0], pf.Rectangle((ox, oy), (ox + length, oy + d)))
        return (ox - 1, oy - 1, ox + length + 1, oy + d + 1)
    if kind == "space":
        w = max(1.0, 2 * limit)
        c.add(layers[0], pf.Rectangle((ox, oy), (ox + length, oy + w)))
        c.add(layers[0], pf.Rectangle((ox, oy + w + d), (ox + length, oy + 2 * w + d)))
        return (ox - 1, oy - 1, ox + length + 1, oy + 2 * w + d + 1)
    inner, outer = layers
    s = max(1.0, 1.2 * limit)
    c.add(inner, pf.Rectangle((ox + d, oy + d), (ox + d + s, oy + d + s)))
    c.add(outer, pf.Rectangle((ox, oy), (ox + s + 2 * d, oy + s + 2 * d)))
    return (ox - 1, oy - 1, ox + s + 2 * d + 1, oy + s + 2 * d + 1)


def _points(marker):
    if hasattr(marker, "hull"):
        return marker.hull
    if hasattr(marker, "edge_pair"):
        return [p for edge in marker.edge_pair for p in edge]
    if hasattr(marker, "edge"):
        return list(marker.edge)
    return []


def test_every_width_space_and_enclosure_rule_fires_exactly_at_its_limit(technology):
    cases = _cases()
    assert len(cases) == 28, f"{len(cases)} width, space and enclosure rules parsed from the deck"
    c = pf.Component("threshold", technology)
    sites = []
    for i, (rule, kind, layers, limit) in enumerate(cases):
        for j, (d, flagged) in enumerate(((round(limit - 0.001, 3), True), (limit, False))):
            box = _site(c, kind, layers, d, j * SITE, -i * SITE, limit)
            sites.append((rule, d, flagged, box))
    c.add("FLOORPLAN", pf.Rectangle((-50, -len(cases) * SITE), (2 * SITE + 50, 50)))
    results = check(c)

    wrong = []
    for rule, d, flagged, (x0, y0, x1, y1) in sites:
        violation = results.violations_by_category.get(rule)
        markers = violation.markers if violation else ()
        hit = any(
            x0 <= x <= x1 and y0 <= y <= y1 for m in markers for x, y in _points(m)
        )
        if hit != flagged:
            wrong.append(f"{rule} at {d} um: {'missed' if flagged else 'false alarm'}")
    assert not wrong, wrong


# ---------------------------------------------------------------------------------------
# The library itself
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", siphopdk.component_names)
def test_every_cell_is_drc_clean_with_its_defaults(name, technology):
    cell = getattr(siphopdk.component, name)(technology=technology)
    assert drc_counts(cell) == {}
