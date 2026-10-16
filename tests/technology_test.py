import matplotlib
import photonforge as pf

import siphopdk

matplotlib.use("Agg")


def test_export(tmp_path):
    tech = siphopdk.siphopdk()
    tech_file = tmp_path / "tech.phf"
    pf.write_phf(tech_file, tech)
    tech_loaded = pf.load_phf(tech_file)["technologies"][0]
    assert tech_loaded.name == tech.name
    assert tech_loaded.version == tech.version
    assert tech_loaded.layers == tech.layers
    assert tech_loaded.ports == tech.ports
    assert tech_loaded.extrusion_specs == tech.extrusion_specs
    assert tech_loaded.background_medium == tech.background_medium
    assert tech_loaded == tech


def test_demo():
    _ = siphopdk.plot_cross_section()


def test_layer_rules_repeat_the_drc_deck():
    """min_width and min_spacing are what PhotonForge routers draw to, so they must equal the deck."""
    deck = siphopdk.DRC_RULES.read_text()
    tech = siphopdk.siphopdk()
    for name, rule_w, rule_s in (("SI", "SI.W.1", "SI.S.1"), ("M1", "M1.W.1", "M1.S.1"),
                                 ("M2", "M2.W.1", "M2.S.1"), ("HTR", "HTR.W.1", "HTR.S.1")):
        spec = tech.layers[name]
        for value, rule in ((spec.min_width, rule_w), (spec.min_spacing, rule_s)):
            line = next(ln for ln in deck.splitlines() if f'output("{rule}"' in ln)
            assert f"({value:.2f}" in line or f"({value:.3f}" in line or f"({value:g}," in line, (
                f"{name}: {value} is not the limit of {rule}: {line}"
            )


def test_metals_conduct_and_are_connected():
    tech = siphopdk.siphopdk()
    conductive = {name for name, spec in tech.layers.items() if spec.conductive}
    assert conductive == {"NPP", "PPP", "VIA1", "VIA2", "VIA3", "M1", "M2"}
    assert not tech.layers["HTR"].conductive, "the heater is a resistor, not a wire"
    layer = {name: tuple(spec.layer) for name, spec in tech.layers.items()}
    pairs = {frozenset(map(tuple, pair)) for pair in tech.connected_layers}
    for a, b in (("VIA1", "M1"), ("VIA2", "M1"), ("M1", "VIA3"), ("VIA3", "M2")):
        assert frozenset((layer[a], layer[b])) in pairs
