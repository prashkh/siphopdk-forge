import warnings

import photonforge as pf
import pytest

import siphopdk


def test_components():
    technology = siphopdk.siphopdk()
    for name in siphopdk.component_names:
        _ = getattr(siphopdk.component, name)(technology=technology)


def test_component_names_are_the_registered_pcells():
    """component_names lists every PCell of the package, and nothing else."""
    import photonforge.parametric_utils as pu

    registered = {key.rsplit(".", 1)[-1] for key in pu._component_registry if key.startswith("siphopdk.")}
    assert set(siphopdk.component_names) == registered
    assert len(siphopdk.component_names) == 27


@pytest.mark.parametrize("name", siphopdk.component_names)
def test_defaults_build_without_warnings(name):
    """A cell built with the default technology unset must not change it or warn."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cell = getattr(siphopdk.component, name)(technology=siphopdk.siphopdk())
    assert cell.technology.name == "SiPhoPDK"
    assert cell.name, f"{name} has no name"


def test_cells_use_the_default_technology_when_it_is_siphopdk():
    tech = siphopdk.siphopdk(sidewall_angle=0.0)
    pf.config.default_technology = tech
    cell = siphopdk.component.straight()
    assert cell.technology is tech or cell.technology == tech


def test_names_are_unique_for_the_defaults():
    technology = siphopdk.siphopdk()
    names = [getattr(siphopdk.component, n)(technology=technology).name for n in siphopdk.component_names]
    assert len(set(names)) == len(names), sorted(names)


def test_photonforge_gui_finds_the_pdk():
    """The component library of the PhotonForge GUI lists the PDK and builds all its cells."""
    from photonforge._live_viewer import pdks

    found = [pdk for pdk in pdks.installed() if pdk["module"] == "siphopdk"]
    if not found:
        pytest.skip("siphopdk-forge is not installed in this environment (running from a checkout)")
    loaded = pdks.load(found[0])
    assert loaded.get("technology") == "SiPhoPDK", loaded.get("error")
    assert not loaded["failed"], loaded["failed"]
    assert len(loaded["entries"]) == len(siphopdk.component_names)
