"""SiPhoPDK rule files, and DRC and LVS with them.

The package ships two KLayout rule files:

* ``DRC_RULES``, ``rules/siphopdk.drc``: the design rule deck, 43 rules.
* ``LVS_RULES``, ``rules/siphopdk.lvs``: the layout-versus-schematic runset.

``drc_check`` and ``lvs_check`` run them through ``photonforge.drc_check`` and
``photonforge.lvs_check``, which need the KLayout application (0.29.12 or newer). The rule
files can also be given to KLayout directly, or chosen in the DRC panel of the PhotonForge GUI.
"""

from pathlib import Path

import photonforge as pf

RULES_DIR = Path(__file__).resolve().parent / "rules"

#: KLayout DRC runset for SiPhoPDK.
DRC_RULES = RULES_DIR / "siphopdk.drc"

#: KLayout LVS runset for SiPhoPDK.
LVS_RULES = RULES_DIR / "siphopdk.lvs"


def drc_check(source, **kwargs):
    """Run the SiPhoPDK design rule deck on a component or a layout file.

    Args:
        source: Component to check, or the path of a GDSII or OASIS file.
        **kwargs: Keyword arguments for ``photonforge.drc_check``, such as ``results_file``
          or ``show_progress``.

    Returns:
        photonforge.verification.DRCResults: Markers grouped by rule. ``is_clean`` is the
        verdict.
    """
    kwargs.setdefault("engine", "klayout")
    return pf.drc_check(source, DRC_RULES, **kwargs)


_LABEL_LAYER = (0, 0)
_PIN_LAYER = (101, 0)
_PIN_E_LAYER = (102, 0)


def _is_route(component):
    return (component.parametric_function or "").startswith("photonforge.parametric.route")


def lvs_pins(design):
    """Pre-export function that names the ports and terminals of every circuit in ``design``.

    Every SiPhoPDK cell draws its own pin boxes and port names. A circuit, a cell that places
    other cells, has none of its own, so the ports it exposes have no name in the layout and
    the runset cannot make pins of them. This function gives each circuit in ``design`` a pin
    box and a text label at each of its optical ports and electrical terminals, on the layers
    the runset reads. Routes from ``photonforge.parametric`` get nothing: the schematic treats
    them as wires.

    Args:
        design: Component to be checked. Its hierarchy decides which cells are circuits.

    Returns:
        Function for ``pre_export_function`` of ``photonforge.lvs_check``.
    """
    from ._recognition import PIN_SIZE, has_pin_marker

    # The export marks a copy of each cell, under its own name, whose references are the
    # original children: a cell is known by its name and its children. The table keeps the
    # children, so that their ids stay valid while the export runs.
    circuits = {}

    def classify(component):
        if _is_route(component):
            return
        children = [reference.component for reference in component.references]
        key = (component.name, tuple(id(child) for child in children))
        if not children or key in circuits:
            return
        circuits[key] = children
        for child in children:
            classify(child)

    classify(design)

    def mark(cell):
        key = (cell.name, tuple(id(reference.component) for reference in cell.references))
        if key not in circuits:
            return
        named = {
            (label.text, round(label.origin[0], 4), round(label.origin[1], 4))
            for labels in cell.labels.values()
            for label in labels
        }
        pins = [(name, port.center, _PIN_LAYER) for name, port in cell.ports.items()
                if has_pin_marker(port)]
        pins += [(name, terminal.center(), _PIN_E_LAYER) for name, terminal in cell.terminals.items()]
        for name, (x, y), layer in pins:
            if (name, round(float(x), 4), round(float(y), 4)) in named:
                continue
            cell.add(layer, pf.Rectangle(center=(x, y), size=(PIN_SIZE, PIN_SIZE)))
            cell.add(_LABEL_LAYER, pf.Label(name, (x, y)))

    return mark


def lvs_check(source, **kwargs):
    """Compare the layout of a component with its own netlist, using the SiPhoPDK runset.

    PhotonForge writes the schematic from the component's connections, and the runset extracts
    the optical and electrical nets from the layout. Every SiPhoPDK cell draws the markers the
    runset reads (a DEVREC outline, a pin box at each port and terminal, and the port names).
    When ``source`` is a component, ``lvs_pins`` adds the same markers to the ports and
    terminals of each circuit in it, so a circuit made of SiPhoPDK cells and PhotonForge routes
    needs no extra preparation.

    Args:
        source: Component to check, or the path of a GDSII or OASIS file (then ``schematic``
          is required).
        **kwargs: Keyword arguments for ``photonforge.lvs_check``, such as ``schematic``,
          ``results_file``, ``pre_export_function`` or ``show_progress``.

    Returns:
        photonforge.verification.LVSResults: ``verdict`` is ``"match"`` when the layout
        connects what the netlist says.
    """
    kwargs.setdefault("engine", "klayout")
    if isinstance(source, pf.Component) and kwargs.get("pre_export_function") is None:
        kwargs["pre_export_function"] = lvs_pins(source)
    return pf.lvs_check(source, LVS_RULES, **kwargs)
