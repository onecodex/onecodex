import click
import re
import os
import shutil
import tempfile
from collections import defaultdict


def auto_detect_pairs(files, prompt):
    """Group paired-end files into tuples in the files list.

    Returns the files list with paired-end files represented as tuples on that list.
    If `prompt` is set to True, the user is asked whether this should happen first.
    """

    # "intelligently" find paired files and tuple them
    paired_files = []
    single_files = set(files)

    for filename in files:
        # convert "read 1" filenames into "read 2" and check that they exist; if they do
        # upload the files as a pair, autointerleaving them
        paired_filename = re.sub(r"([._][Rr])1([._][\w.]+)$", r"\g<1>2\g<2>", filename)
        paired_filename = re.sub(r"([._])1([._][\D._]+)$", r"\g<1>2\g<2>", paired_filename)

        # we don't necessary need the R2 to have been passed in; we infer it anyways
        if paired_filename != filename and os.path.exists(paired_filename):
            if not prompt and paired_filename not in single_files:
                # if we're not prompting, don't automatically pull in files
                # not in the list the user passed in
                continue

            paired_files.append((filename, paired_filename))

            if paired_filename in single_files:
                single_files.remove(paired_filename)

            single_files.remove(filename)

    auto_pair = True

    if prompt and len(paired_files) > 0:
        pair_list = ""
        for p in paired_files:
            pair_list += "\n  {}  &  {}".format(os.path.basename(p[0]), os.path.basename(p[1]))

        answer = click.confirm(
            "It appears there are {n_paired_files} paired files (of {n_files} total):{pair_list}\nInterleave them after upload?".format(
                n_paired_files=len(paired_files) * 2,
                n_files=len(paired_files) * 2 + len(single_files),
                pair_list=pair_list,
            ),
            default="Y",
        )

        if not answer:
            auto_pair = False

    if auto_pair:
        return paired_files + list(single_files)
    else:
        return files


def _find_multilane_groups(files):
    """Find a list of multilane file groups eligible for concatenation.

    The files are grouped based on filename (e.g. `Sample_R1_L001.fq`, `Sample_R1_L002.fq`).
    If there is a gap in the sequence (e.g. [`Sample_R1_L001.fq`, `Sample_R1_L003.fq`]), the group
    is skipped. If there is a mismatch in forward and reverse file sequence (e.g.
    [(`Sample_R1_L001.fq`, `Sample_R2_L001.fq`), `Sample_R2_L002.fq`]), the group is skipped.
    If the sequence doesn't begin with `L001`, the group is skipped.

    This function assumes that the paired-end file tuples on the list are properly matched.

    The result is a list of lists, with each nested list representing a single multilane
    file group concisting of either string filenames (for single read files) or tuples
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


def concatenate_multilane_files(files, prompt):
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
        target_path = os.path.join(tempfile.gettempdir(), target_file_name)

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
