import { loadPyodide } from "pyodide";

async function main() {
  let pyodide = await loadPyodide();
  pyodide.FS.mkdir("/source_code");
  pyodide.FS.mount(pyodide.FS.filesystems.NODEFS, { root: "../../" }, "/source_code");

  await pyodide.loadPackage("micropip");

  // Installing libs from pyproject.toml
  //
  // Some are just DQL libs (like ruff) and are not needed for running tests.
  //
  // There is a special hack for scikit-bio. It has an int overflow bug and the fix
  // is not released yet.
  await pyodide.runPythonAsync(`
    import tomllib
    from pathlib import Path

    import micropip

    root_path = Path("/source_code")
    project_path = root_path / "pyproject.toml"
    with open(project_path, "rb") as h:
        project_cfg = tomllib.load(h)

    PACKAGES_TO_IGNORE = {
        "coverage",
        "pytest-cov",
        "pdfplumber",
        "ruff",
        "scikit-bio",
    }

    def should_install(package):
        for x in PACKAGES_TO_IGNORE:
            if package.startswith(x):
                return False
        return True

    await micropip.install(
        "https://static.onecodex.com/app/pyodide/v0.29.0/"
        "scikit_bio-0.7.1.post7-cp313-cp313-pyodide_2025_0_wasm32.whl",
        deps=False,
    )

    for package in project_cfg["dependency-groups"]["dev"]:
        if should_install(package):
            await micropip.install(package)

    for package in project_cfg["project"]["dependencies"]:
        if should_install(package):
            await micropip.install(package)

    for package in project_cfg["project"]["optional-dependencies"]["all"]:
        if should_install(package):
            await micropip.install(package)
  `);

  // Mocking some libraries
  pyodide.runPython(`
    import sys, types

    class NoOp:
        def __getattr__(self, name):
            return self

        def __setattr__(self, *args, **kwargs):
            pass

    class MockModule(types.ModuleType):
        def __getattr__(self, name):
            return NoOp()

    mocked_modules = [
        "array_api_compat",
        "biom",
        "filelock",
        "h5py",
        "nbformat",
        "sentry_sdk",
        "sqlite3",
    ]

    for module in mocked_modules:
        sys.modules[module] = MockModule(module)

    sys.modules["biom"].Table = NoOp
  `);

  // Actually running the tests
  await pyodide.runPythonAsync(`
    import sys
    import os
    from pathlib import Path

    sys.path.append("/source_code")
    os.chdir("/source_code/")

    Path("__init__.py").touch()

    import pytest
    exit_code = pytest.main([
        "-x",
        "/source_code/tests/test_dataframes.py",
        "/source_code/tests/test_collection.py",
        "/source_code/tests/test_custom_plots.py",
        "/source_code/tests/test_viz.py",
        "/source_code/tests/test_stats.py",
    ])

    assert exit_code == 0
  `)
}

main();
