import { loadPyodide } from "pyodide";

async function main() {
  let pyodide = await loadPyodide();
  pyodide.FS.mkdir("/source_code");
  pyodide.FS.mount(pyodide.FS.filesystems.NODEFS, { root: "../" }, "/source_code");

  await pyodide.loadPackage("micropip");
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

    await micropip.install("h5py")
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
        "filelock",
        "biom",
        "sentry_sdk",
        "sqlite3",
        "array_api_compat",
        "nbformat",
    ]

    for module in mocked_modules:
        sys.modules[module] = MockModule(module)

    sys.modules["biom"].Table = NoOp
  `);

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
        "/source_code/tests/test_analyses.py",
        "/source_code/tests/test_dataframes.py",
        "/source_code/tests/test_collection.py",
        "/source_code/tests/test_custom_plots.py",
        "/source_code/tests/test_viz.py",
        "/source_code/tests/test_distance.py",
        "/source_code/tests/test_stats.py",
    ])

    assert exit_code == 0
  `)
}

main();
