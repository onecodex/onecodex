import click
import re
import sys
import os
import shutil
import logging
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

# Captures the ordinal and the parts before and after it
# (. or _ followed by num followed by . or _ and non-digits)
ORDINAL_REV_PATTERN = r"(?P<pre>[._])(?P<ordinal>\d)(?P<post>[._][\D._]+)$"
ORDINAL_MULTI_REV_PATTERN = r"(?P<pre>[._])(?P<ordinal>\d+)(?P<post>[._][\D._]+)$"

# Captures parts before and after paired ordinal (as above, but includes R1, r1, R2, r2)
PAIRED_ORDINAL_REV_PATTERN = r"(?P<pre>[._][Rr])(?P<ordinal>\d)(?P<post>[._][\w.]+)$"

# Captures the sequencing lane number
LANE_PATTERN = re.compile(r"[._]L(?P<lane>\d+)(?=[._])")

log = logging.getLogger("onecodex")


@dataclass(frozen=True)
class PlannedSample:
    """One sample to upload, and the files it will be assembled from.

    `forward` holds the files that make up the sample, in order; more than one means they
    are concatenated. `reverse` holds the corresponding files for the second read of a
    paired end sample, and is None for a single ended one.
    """

    forward: tuple[str, ...]
    reverse: tuple[str, ...] | None = None

    @property
    def files(self) -> tuple[str, ...]:
        return self.forward + (self.reverse or ())

    @property
    def is_paired(self) -> bool:
        return self.reverse is not None

    @property
    def is_concatenated(self) -> bool:
        return len(self.forward) > 1


def _replace_filename_ordinal(filename, replacement, multi_digit=False):
    """Replace file[_.]num[._] with file[_.]replacement[._].

    If `multi_digit` is set to True, num may be a multi digit number
    """
    replace_pattern = rf"\g<pre>{replacement}\g<post>"
    regex = ORDINAL_MULTI_REV_PATTERN if multi_digit else ORDINAL_REV_PATTERN
    return re.sub(regex, replace_pattern, filename)


def _replace_paired_filename_ordinal(filename, replacement):
    """Replace file[_.][Rr]?num[._] with file[_.][Rr]?replacement[._]."""
    first_pass = _replace_filename_ordinal(filename, replacement)
    replace_pattern = rf"\g<pre>{replacement}\g<post>"
    return re.sub(PAIRED_ORDINAL_REV_PATTERN, replace_pattern, first_pass)


def _lane_number(path: str) -> int | None:
    """Return the sequencing lane a file belongs to, if its name carries one."""
    match = LANE_PATTERN.search(path)
    return int(match.group("lane")) if match else None


def _assembled_name(path: str) -> str:
    """Return the filename a sample gets once its parts have been joined together."""
    name = LANE_PATTERN.sub("", path)
    return re.sub(ORDINAL_MULTI_REV_PATTERN, r"\g<post>", name)


def _ont_run_on_disk(filename: str) -> list[str]:
    """Return the run of ONT files on disk starting at ordinal 0, stopping at the first gap."""
    run = []
    idx = 0
    while True:
        sibling = _replace_filename_ordinal(filename, idx, multi_digit=True)
        if not os.path.isfile(sibling):
            return run
        run.append(sibling)
        idx += 1


def _ont_sequence(group: set[str], prompt: bool) -> list[str] | None:
    """Return a group's chunks in order, or None if they are not a whole sample.

    A sample's chunks are a contiguous run starting at ordinal 0. When prompting, the run
    may be completed from disk; otherwise it has to be complete among the files passed in.
    """
    any_file = next(iter(group))
    if prompt:
        sequence = _ont_run_on_disk(any_file)
    else:
        sequence = [
            _replace_filename_ordinal(any_file, idx, multi_digit=True) for idx in range(len(group))
        ]

    # a single chunk is a whole file already, so there is nothing to assemble
    if len(sequence) < 2 or not group <= set(sequence):
        return None
    return sequence


def _plan_ont_samples(files: Sequence[str], prompt: bool) -> tuple[list[PlannedSample], list[str]]:
    """Claim the files that make up ONT samples split across numbered chunks.

    A sample's chunks are a contiguous run starting at ordinal 0. The whole run has to be
    accounted for: files that cannot be part of one are left for the other planners, so
    that a pair of Illumina reads is never mistaken for a partial run.
    """
    candidates = defaultdict(set)
    for filename in files:
        # without an ordinal the substitution is a no-op and the file matches itself
        if not re.search(ORDINAL_MULTI_REV_PATTERN, filename):
            continue

        # the directory is part of the key: two samples may share a filename
        base = re.sub(ORDINAL_MULTI_REV_PATTERN, r"\g<post>", filename)
        candidates[base].add(filename)
        if prompt:
            # only pull in files the user did not name if we can ask them about it
            run = _ont_run_on_disk(filename)
            if filename in run:
                candidates[base].update(run)

    samples, claimed = [], set()
    for group in candidates.values():
        sequence = _ont_sequence(group, prompt)
        if sequence is None:
            continue
        samples.append(PlannedSample(forward=tuple(sequence)))
        claimed.update(group)

    return samples, [f for f in files if f not in claimed]


def _plan_paired_samples(
    files: Sequence[str], prompt: bool
) -> tuple[list[PlannedSample], list[str]]:
    """Claim the files that make up paired end samples.

    The mate of a file is found by substituting the read number into its name, so the two
    always live in the same directory. When prompting, a mate that was not passed in may
    still be picked up from disk.
    """
    named = set(files)
    samples, claimed = [], set()

    for filename in files:
        if filename in claimed:
            continue

        r1 = _replace_paired_filename_ordinal(filename, "1")
        r2 = _replace_paired_filename_ordinal(filename, "2")

        # a file with any other ordinal substitutes down to the same two names without
        # being either of them
        if r1 == r2 or filename not in (r1, r2):
            continue
        if not (os.path.isfile(r1) and os.path.isfile(r2)):
            continue

        mate = r2 if filename == r1 else r1
        if not prompt and mate not in named:
            continue

        samples.append(PlannedSample(forward=(r1,), reverse=(r2,)))
        claimed.update((r1, r2))

    return samples, [f for f in files if f not in claimed]


def _merge_lanes(samples: list[PlannedSample]) -> list[PlannedSample]:
    """Combine samples that are the same library sequenced across several lanes.

    The lanes have to run from 1 without a gap, and every lane has to agree about whether
    the sample is paired end, otherwise the group is left alone.
    """
    groups = defaultdict(list)
    for sample in samples:
        if _lane_number(sample.forward[0]) is None:
            continue
        groups[LANE_PATTERN.sub("", sample.forward[0])].append(sample)

    merged, consumed = [], set()
    for group in groups.values():
        if len(group) < 2:
            continue

        group.sort(key=lambda s: _lane_number(s.forward[0]))
        lanes = [_lane_number(s.forward[0]) for s in group]
        if lanes != list(range(1, len(group) + 1)):
            continue
        if len({s.is_paired for s in group}) > 1:
            continue

        forward = tuple(f for s in group for f in s.forward)
        reverse = tuple(f for s in group for f in s.reverse) if group[0].is_paired else None
        merged.append(PlannedSample(forward=forward, reverse=reverse))
        consumed.update(group)

    return merged + [s for s in samples if s not in consumed]


def plan_uploads(files: Sequence[str], prompt: bool) -> list[PlannedSample]:
    """Work out which samples the given files make up.

    Each file belongs to exactly one sample. ONT chunks are claimed first, so that two
    chunks of one sample are never mistaken for a pair of Illumina reads; what is left is
    paired up, and anything still unclaimed is a sample on its own. Finally, samples that
    are the same library across several lanes are combined.
    """
    ont_samples, rest = _plan_ont_samples(files, prompt)
    paired_samples, rest = _plan_paired_samples(rest, prompt)
    singles = [PlannedSample(forward=(filename,)) for filename in rest]

    return _merge_lanes(ont_samples + paired_samples + singles)


def _describe(sample: PlannedSample) -> str:
    """Return a short phrase saying what will be done to a sample's files."""
    steps = []
    if sample.is_concatenated:
        if _lane_number(sample.forward[0]) is None:
            steps.append(f"concatenate {len(sample.forward)} files")
        else:
            steps.append(f"concatenate {len(sample.forward)} lanes")
    if sample.is_paired:
        steps.append("interleave")

    return ", then ".join(steps) if steps else "upload as-is"


def describe_plan(samples: Sequence[PlannedSample], named: set[str]) -> None:
    """Print what each sample will be built from, marking files found on disk."""

    def _label(path):
        return f"{path}" if path in named else f"{path} *"

    click.echo(
        click.wrap_text(
            "One Codex stores each sample as a single file, so files belonging to the same "
            "sample are joined together before they are uploaded: paired end reads are "
            "interleaved, and files split into numbered chunks or across sequencing lanes "
            "are concatenated. Which files belong together is worked out from their names, "
            "so please check this over before continuing.",
            width=min(shutil.get_terminal_size().columns, 88),
        )
    )
    click.echo()
    click.echo("Planned uploads:\n")

    for position, sample in enumerate(samples, start=1):
        click.echo(f"{position:>2}. {_describe(sample)}")
        if sample.is_paired:
            rows = [
                " + ".join(_label(f) for f in sample.forward),
                " + ".join(_label(f) for f in sample.reverse),
            ]
        else:
            rows = [_label(f) for f in sample.forward]
        for n, row in enumerate(rows, start=1):
            prefix = "└──" if n == len(rows) else "├──"
            click.echo(f"    {prefix} {row}")
        click.echo()

    n_files = sum(len(s.files) for s in samples)
    if any(f not in named for s in samples for f in s.files):
        click.echo("* not specified on the command line; found alongside the files that were\n")
    click.echo(f"{len(samples)} sample(s) from {n_files} file(s).\n")


def confirm_plan(samples: Sequence[PlannedSample], named: set[str]) -> bool:
    """Show the plan and ask whether to go ahead with it.

    Returns False if the files should be uploaded exactly as they were given instead.
    """
    describe_plan(samples, named)

    answer = click.prompt(
        "[Y]es, upload as planned; [n]o, upload each file I specified as a separate sample;"
        " [c]ancel",
        type=click.Choice(["Y", "n", "c"], case_sensitive=False),
        default="Y",
    )

    if answer.lower() == "c":
        click.echo("Upload canceled")
        sys.exit(0)

    return answer.lower() == "y"


def _concatenate(paths: Sequence[str], tempdir: str, position: int, suffix: str) -> str:
    """Join files together into `tempdir`, keeping the name the sample will upload under.

    Each sample writes into its own subdirectory, so two samples whose files share a name
    but live in different directories do not overwrite one another.
    """
    sample_dir = os.path.join(tempdir, f"{position}{suffix}")
    os.makedirs(sample_dir, exist_ok=True)
    target = os.path.join(sample_dir, os.path.basename(_assembled_name(paths[0])))

    log.info(f"Concatenating to {target}")
    with open(target, "wb") as outf:
        for path in paths:
            with open(path, "rb") as inf:
                shutil.copyfileobj(inf, outf)
    return target


def materialize_plan(samples: Sequence[PlannedSample], tempdir: str) -> list[str | tuple[str, str]]:
    """Build the files the plan describes, and return them ready to upload."""
    uploads = []
    for position, sample in enumerate(samples, start=1):
        if sample.is_concatenated:
            forward = _concatenate(sample.forward, tempdir, position, "")
            reverse = (
                _concatenate(sample.reverse, tempdir, position, "r") if sample.is_paired else None
            )
        else:
            forward = sample.forward[0]
            reverse = sample.reverse[0] if sample.is_paired else None

        uploads.append((forward, reverse) if reverse else forward)
    return uploads
