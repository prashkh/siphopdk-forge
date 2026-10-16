"""Publish SiPhoPDK as a public PDA library, so that it shows in the PhotonForge cloud GUI.

    python make_pda_library.py

The library version is the technology version. A version that already exists is not published
again.
"""

import shutil

import photonforge as pf
from photonforge import pda

import siphopdk


def ignore_cache(path, names):
    return ["__pycache__"] if "__pycache__" in names else []


def create_library():
    name = "SiPhoPDK"
    version = siphopdk.siphopdk().version
    description = "SiPhoPDK: synthetic 220 nm SOI + SiN O-band demonstration PDK"

    for lib in pda.list_libraries(name):
        if lib["version"] == version:
            print("Library already exists: " + str(lib))
            return

    print(f"Creating library {name!r} - version {version!r}", flush=True)
    project = pda.create_project(
        name=name,
        description=description,
        visibility="public",
        role="viewer",
        create_template=False,
    )

    # Add sources
    shutil.copytree(
        "./siphopdk",
        project.module_path / project.module_name,
        dirs_exist_ok=True,
        ignore=ignore_cache,
    )

    project.save_module()
    module = project.import_module(None)[project.module_name]

    tech = module.siphopdk()
    pf.config.default_technology = tech

    # Pass the technology, so that it is set on updates
    components = [
        getattr(module.component, n)(technology=tech) for n in siphopdk.component_names
    ]

    update_config = True
    for component in components:
        print("Adding", repr(component.name), flush=True)
        project.add(component, update_existing_dependencies=False, update_config=update_config)
        update_config = False

    project.add_version(version)

    print(f"Done: {project!r}")


if __name__ == "__main__":
    pda.init()

    try:
        create_library()
    finally:
        pda.stop()
