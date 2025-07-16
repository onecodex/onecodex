# Test import speed of the module. This is important as it affects click auto-complete
# and our general CLI startup time.
import os
import subprocess
import time

import pytest

from onecodex import Cli


@pytest.mark.parametrize(
    "import_command,package_max_import_times,total_time_secs",
    [
        # Fast imports
        ("import onecodex", {"onecodex": 0.30}, 0.30),
        ("from onecodex import Cli", {"onecodex": 0.30, "onecodex.cli": 0.20}, 0.30),
        ("from onecodex import Api", {"onecodex": 0.30, "onecodex.api": 0.20}, 0.30),
        # Full startup of the CLI (prints help message)
        ('import onecodex; onecodex.cli.onecodex(["--help"])', {"onecodex": 0.30}, 0.30),
        # Our slow imports should be lazy and still import fast
        (
            "from onecodex.viz import VizPCAMixin",
            {
                "onecodex": 0.30,
                "onecodex.viz": 0.6,
                "onecodex.viz._pca": 0.01,
                "onecodex.viz._distance": 0.01,
            },
            0.6,
        ),
        (
            "from onecodex.analyses import AnalysisMixin",
            {"onecodex": 0.30, "onecodex.analyses": 0.30},
            0.30,
        ),
    ],
)
def test_import_speed(import_command, package_max_import_times, total_time_secs):
    """Test the import speed for key onecodex modules

    Note we test the `package`, so `onecodex` will test *everything*
    while `onecodex.cli` will test the incremental impact of loads the cli.py
    module.
    """
    from pathlib import Path

    try:
        import pandas as pd  # noqa
    except ImportError:
        assert False, "Must run test in an environment with onecodex[all] dependencies"

    # Get the virtualenv and try to activate it
    activate_path = Path(os.environ["VIRTUAL_ENV"]).joinpath("bin/activate")
    cmd = [
        "/bin/bash",
        "-c",
        "source {} && python -X importtime -c '{}'".format(activate_path, import_command),
    ]
    cmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmd_out, cmd_err = cmd.communicate()

    # Format is like `import time:       227 |     112546 | onecodex`
    # Convert to a list of lists and skip anything with a header value (`self [us]`)
    table = [
        [y.strip() for y in x.replace("import time:", "").split("|")]
        for x in cmd_err.decode().split("\n")
        if x and "self [us]" not in x
    ]
    df = pd.DataFrame(table)
    df.columns = ["time", "cumulative", "package"]
    df.time = df.time.astype(int)
    df.cumulative = df.cumulative.astype(int)

    # Check times
    for package, max_time_secs in package_max_import_times.items():
        cumulative_import_time_us = df[df.package == package].cumulative.values[0]
        cumulative_import_time_secs = cumulative_import_time_us / 1e6
        assert cumulative_import_time_secs < max_time_secs

    # Check max time
    assert max(df.cumulative) / 1e6 < total_time_secs


def test_cli_speed(runner, api_data, mocked_creds_file):
    """Test that loading the CLI is fast, i.e., it doesn't load all the extensions"""
    start_time = time.time()
    result = runner.invoke(Cli, ["analyses"])
    elapsed = time.time() - start_time
    assert result.exit_code == 0
    assert elapsed <= 0.2
