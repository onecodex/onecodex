import click
import re
import sys
import os
import shutil
import logging
from collections import defaultdict
from collections.abc import Sequence

# Captures the ordinal and the parts before and after it
# (. or _ followed by num followed by . or _ and non-digits)
ORDINAL_REV_PATTERN = r"(?P<pre>[._])(?P<ordinal>\d)(?P<post>[._][\D._]+)$"
ORDINAL_MULTI_REV_PATTERN = r"(?P<pre>[._])(?P<ordinal>\d+)(?P<post>[._][\D._]+)$"

# Captures parts before and after paired ordinal (as above, but includes R1, r1, R2, r2)
PAIRED_ORDINAL_REV_PATTERN = r"([._][Rr])\d([._][\w.]+)$"

log = logging.getLogger("onecodex")


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
    replace_pattern = rf"\g<1>{replacement}\g<2>"
    return re.sub(PAIRED_ORDINAL_REV_PATTERN, replace_pattern, first_pass)


def _ont_sequence_on_disk(filename: str) -> list[str]:
    """Return the run of ONT files on disk starting at ordinal 0, stopping at the first gap.

    Returns an empty list if the run does not reach `filename`, which is then not part of a
    complete sequence: the files before the gap belong to no sample we could assemble.
    """
    sequence = []
    idx = 0
    while True:
        sibling = _replace_filename_ordinal(filename, idx, multi_digit=True)
        if not os.path.isfile(sibling):
            break
        sequence.append(sibling)
        idx += 1

    return sequence if filename in sequence else []


def prompt_user_for_concatenation(ont_groups: dict[str, set[str]], passed_files: set[str]) -> bool:
    """Prompt user to determine whether ONT files should be concatenated.

    Files in `ont_groups` that are not in `passed_files` were found on disk rather than
    given on the command line, and are marked as such.
    """

    n_files = sum([len(x) for x in ont_groups.values()])
    n_unspecified = sum(1 for files in ont_groups.values() for f in files if f not in passed_files)

    message = f"It appears there are {len(ont_groups)} sample(s) split across {n_files} individual file(s)."
    if n_unspecified:
        message += (
            f"\n{n_unspecified} of them were not specified on the command line; "
            "they were found alongside the files that were."
        )

    answer = click.prompt(
        message + "\n\nWould you like to concatenate files by sample?"
        "\n\n[Y]es (recommended); [n]o, upload each file I specified as a separate sample;"
        " [d]isplay files; [c]ancel",
        type=click.Choice(["Y", "n", "d", "c"]),
        default="Y",
    )

    if answer[0] == "Y":
        return True
    elif answer[0] == "n":
        return False
    elif answer[0] == "c":
        click.echo("Upload canceled")
        sys.exit(0)
    elif answer[0] == "d":

        def _ordinal(filename: str) -> int:
            # sorting the filenames directly would put sample_10 before sample_2
            return int(re.search(ORDINAL_MULTI_REV_PATTERN, filename).group("ordinal"))

        for group_name, group_files in ont_groups.items():
            click.echo(f"{os.path.basename(group_name)}")
            ordered = sorted(group_files, key=_ordinal)
            for n, group_file in enumerate(ordered, start=1):
                if n == len(ordered):
                    prefix = "└──"
                else:
                    prefix = "├──"
                marker = "" if group_file in passed_files else " *"
                click.echo(f"{prefix} {group_file}{marker}")
            click.echo()
        if n_unspecified:
            click.echo("* not specified on the command line; found alongside the other files")
            click.echo()
        return prompt_user_for_concatenation(ont_groups, passed_files)
    else:
        click.echo(f"Unknown option: {answer}")
        return prompt_user_for_concatenation(ont_groups, passed_files)


def concatenate_ont_groups(
    files: Sequence[str], prompt: bool, tempdir: str
) -> tuple[list[str], list[str]]:
    """Concatenate ONT split files.

    Returns `(concatenated, remaining)`.

    If `prompt` is set to True, the rest of a sample's file sequence may be picked up from
    disk even if it was not passed in. Such files are marked as not specified in the prompt.
    """
    passed_files = set(files)
    single_files = set(files)
    concatenated = []
    ont_groups = defaultdict(set)
    auto_group = True

    for filename in files:
        # without an ordinal, the substitution below is a no-op and the file
        # would be matched against itself
        if not re.search(ORDINAL_MULTI_REV_PATTERN, filename):
            continue

        ont_zero_filename = _replace_filename_ordinal(filename, "0", multi_digit=True)
        if os.path.isfile(ont_zero_filename):
            # strip the ordinal and preceding . or _
            base_filename = re.sub(ORDINAL_MULTI_REV_PATTERN, r"\g<post>", filename)
            base_filename = os.path.join(tempdir, os.path.basename(base_filename))

            ont_groups[base_filename].add(filename)
            if prompt:
                # if we're not prompting, don't automatically pull in files
                # not in the list the user passed in
                ont_groups[base_filename].update(_ont_sequence_on_disk(filename))

    if not ont_groups:
        return concatenated, list(files)

    # a lone file is a whole sample already; there is nothing to concatenate
    if len(ont_groups) == 1 and sum(len(x) for x in ont_groups.values()) == 1:
        auto_group = False
    elif prompt:
        auto_group = prompt_user_for_concatenation(ont_groups, passed_files)
    else:
        auto_group = True

    if not auto_group:
        return concatenated, list(files)

    # Ensure there is no gap in the file sequences
    for base_ont_filename, group_files in ont_groups.items():
        ont_file = next(iter(group_files))
        expected_sequence = [
            _replace_filename_ordinal(ont_file, idx, multi_digit=True)
            for idx in range(len(group_files))
        ]
        full_sequence = True
        for expected_file in expected_sequence:
            if expected_file not in group_files:
                log.warning(
                    "Detected a gap in the ONT file sequence for "
                    f"{os.path.basename(base_ont_filename)}, missing file:"
                    f" {os.path.basename(expected_file)}. Skipping concatenation"
                )
                full_sequence = False
                break

        if not full_sequence:
            continue

        log.info(f"Concatenating to {base_ont_filename}")
        with open(base_ont_filename, "wb") as outf:
            for ont_filename in expected_sequence:
                with open(ont_filename, "rb") as inf:
                    shutil.copyfileobj(inf, outf)
                single_files.discard(ont_filename)
        concatenated.append(base_ont_filename)
    return concatenated, list(single_files)


def auto_detect_illumina_pairs(files: Sequence[str], prompt: bool) -> list[str | tuple[str, str]]:
    """Group paired-end files in the files list.

    Returns the files list with paired-end files represented as tuples on that list.
    If `prompt` is set to True, the user is asked whether this should happen first, and
    the mate of a paired file may be picked up from disk even if it was not passed in.
    Such files are marked as not specified in the prompt.
    """

    passed_files = set(files)

    # files left ungrouped
    single_files = set(files)

    # "intelligently" grouped paired-end files
    pairs = []

    for filename in files:
        if filename not in single_files:
            # filename may have been already removed as a pair
            continue

        paired_r1_filename = _replace_paired_filename_ordinal(filename, "1")
        paired_r2_filename = _replace_paired_filename_ordinal(filename, "2")

        if (
            paired_r1_filename != paired_r2_filename
            # a file with any other ordinal (e.g. `sample_3.fq`) substitutes down to the
            # same two names, but is not itself a mate
            and filename in (paired_r1_filename, paired_r2_filename)
            and os.path.isfile(paired_r1_filename)
            and os.path.isfile(paired_r2_filename)
        ):
            other_paired_file = (
                paired_r2_filename if filename == paired_r1_filename else paired_r1_filename
            )
            # we don't necessary need the other paired to have been passed in; we infer it anyways
            if not prompt and other_paired_file not in single_files:
                # if we're not prompting, don't automatically pull in files
                # not in the list the user passed in
                continue

            pairs.append((paired_r1_filename, paired_r2_filename))
            if paired_r1_filename in single_files:
                single_files.remove(paired_r1_filename)
            if paired_r2_filename in single_files:
                single_files.remove(paired_r2_filename)

    auto_pair = True
    if prompt and pairs:

        def _label(filename):
            # mark files we found on disk but that were not passed on the command line
            marker = "" if filename in passed_files else " *"
            return f"{os.path.basename(filename)}{marker}"

        pair_list = ""
        for pair in pairs:
            pair_list += f"\n  {_label(pair[0])}  &  {_label(pair[1])}"
        if any(f not in passed_files for pair in pairs for f in pair):
            pair_list += "\n* not specified on the command line; found alongside its mate"

        summary = "is 1 pair" if len(pairs) == 1 else f"are {len(pairs)} pairs"
        if len(single_files) == 1:
            summary += ", and 1 other file"
        elif single_files:
            summary += f", and {len(single_files)} other files"

        answer = click.prompt(
            f"It appears there {summary}:{pair_list}"
            "\n\nWould you like to interleave each pair?"
            "\n\n[Y]es (recommended); [n]o, upload each file I specified as a separate sample;"
            " [c]ancel",
            type=click.Choice(["Y", "n", "c"], case_sensitive=False),
            default="Y",
        )

        if answer[0].lower() == "c":
            click.echo("Upload canceled")
            sys.exit(0)
        elif answer[0].lower() == "n":
            auto_pair = False

    if auto_pair:
        return pairs + list(single_files)
    else:
        return list(files)


def _find_multilane_groups(files):
    """Find a list of multilane file groups eligible for concatenation.

    The files are grouped based on filename (e.g. `Sample_R1_L001.fq`, `Sample_R1_L002.fq`).
    If there is a gap in the sequence (e.g. [`Sample_R1_L001.fq`, `Sample_R1_L003.fq`]), the group
    is skipped. If there is a mismatch in forward and reverse file sequence (e.g.
    [(`Sample_R1_L001.fq`, `Sample_R2_L001.fq`), `Sample_R2_L002.fq`]), the group is skipped.
    If the sequence doesn't begin with `L001`, the group is skipped.

    This function assumes that the paired-end file tuples on the list are properly matched.

    The result is a list of lists, with each nested list representing a single multilane
    file group consisting of either string filenames (for single read files) or tuples
    (for paired-end reads). The files are in proper order, concatenation-ready.
    """

    pattern_multilane = re.compile(r"[._]L(\d+)[._]")
    pattern_pair_lane_combo = re.compile(r"([._][rR][12])?[._]L\d+[._]([rR][12])?")

    def _group_for(file_path):
        """Create group names by removing Lx and Rx elements from the filename."""
        return re.sub(pattern_pair_lane_combo, "", os.path.basename(file_path))

    def _create_group_map(elem_list, paired):
        """Create multilane file groups with elements in proper order based on file list."""
        # Create groups for the multilane files
        group_map = defaultdict(list)
        for elem in elem_list:
            search_elem = elem if not paired else elem[0]
            if pattern_multilane.search(search_elem):
                group = _group_for(search_elem)
                group_map[group].append(elem)

        # Only multifile groups are returned
        return {
            group: sorted(elems, key=lambda x: x[0] if paired else x)
            for group, elems in group_map.items()
            if len(elems) > 1
        }

    def _with_gaps_removed(group_map, paired):
        """Return a new map having groups with gaps in elements removed."""
        gapped_groups = set()
        for group, elems in group_map.items():
            # Verify we're getting 1, 2, 3, ...
            expected_sequence = list(range(1, len(elems) + 1))
            if paired:
                fwd_nums = [
                    int(pattern_multilane.search(se).group(1)) for se in [fwd for fwd, _ in elems]
                ]
                rev_nums = [
                    int(pattern_multilane.search(se).group(1)) for se in [rev for _, rev in elems]
                ]
                if fwd_nums != expected_sequence or rev_nums != expected_sequence:
                    gapped_groups.add(group)
            else:
                nums = [int(pattern_multilane.search(se).group(1)) for se in elems]
                if nums != expected_sequence:
                    gapped_groups.add(group)

        return {group: elems for group, elems in group_map.items() if group not in gapped_groups}

    single_files = [f for f in files if isinstance(f, str)]
    paired_files = [f for f in files if isinstance(f, tuple)]

    multilane_pairs = _create_group_map(paired_files, paired=True)
    multilane_singles = _create_group_map(single_files, paired=False)

    # Search for unmatched files for paired end multilane files and remove offending groups,
    # e.g. [(Sample_R1_L001.fq, Sample_R2_L001.fq), Sample_R2_L002.fq]
    for filename in single_files:
        if pattern_multilane.search(filename):
            group = _group_for(filename)
            if group in multilane_pairs:
                del multilane_pairs[group]

    # Remove groups with gaps, e.g. [`Sample_R1_L001.fq`, `Sample_R1_L003.fq`]
    multilane_pairs = _with_gaps_removed(multilane_pairs, paired=True)
    multilane_singles = _with_gaps_removed(multilane_singles, paired=False)

    multilane_groups = list(multilane_singles.values())
    multilane_groups.extend(list(multilane_pairs.values()))

    return multilane_groups


def concatenate_multilane_files(files, prompt, tempdir):
    """Concatenate multilane files before uploading.

    The files are grouped based on filename. If `prompt` is set to True, the user
    is asked whether this should happen first.

    The concatenated files replace the matched sequence files. They're put in a temporary
    directory and overwrite any existing files.

    Returns a new list with multilane groups replaced with a path to the single
    concatenated file.
    """

    def _concatenate_group(group, first_elem):
        """Concatenate all the files on the list and return the target file path."""
        target_file_name = re.sub(pattern_lane_num, r"\1", os.path.basename(first_elem))
        target_path = os.path.join(tempdir, os.path.basename(target_file_name))

        # Overwriting all files by default
        with open(target_path, "wb") as outf:
            for fname in group:
                with open(fname, "rb") as inf:
                    # TODO: check for newline at the end of file first?
                    shutil.copyfileobj(inf, outf)
        return target_path

    groups = _find_multilane_groups(files)

    if not groups:
        return files

    perform_concat = True
    if prompt:
        answer = click.confirm(
            "This data appears to have been split across multiple sequencing lanes.\nConcatenate lanes before upload?",
            default="Y",
        )
        if not answer:
            perform_concat = False

    if not perform_concat:
        return files

    files = files[:]
    pattern_lane_num = re.compile(r"[._]L\d+([._])")

    for group in groups:
        # The groups considered here will already have more than 1 element
        first_elem = group[0]
        if isinstance(first_elem, tuple):
            concat_fwd = _concatenate_group([fwd for fwd, _ in group], first_elem[0])
            concat_rev = _concatenate_group([rev for _, rev in group], first_elem[1])
            files.append((concat_fwd, concat_rev))
        elif isinstance(first_elem, str):
            concat = _concatenate_group(group, first_elem)
            files.append(concat)

        for elem in group:
            files.remove(elem)

    return files
