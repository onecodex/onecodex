#!/usr/bin/env python
"""
No dependency Python script for joining two paired end FASTQ files.

Auto-detects gzip'd files, offers header checking via a --strict flag,
and supports output to STDOUT a gzip'd FASTQ or an uncompressed FASTQ (--uncompressed flag).
"""

import click
import gzip

from typing import TextIO, cast


def enforce_headers(f1_header: str, f2_header: str):
    if f1_header[0] != "@" or f2_header[0] != "@":
        raise Exception("Invalid input FASTQ files.")
    if f1_header.strip().split(" ")[0] != f2_header.strip().split(" ")[0]:
        raise Exception("Input FASTQ files do not share headers. Check and re-run w/o --strict.")


def flush_buffer(
    f_out: TextIO, f1_lines: list[str], f2_lines: list[str]
) -> tuple[list[str], list[str]]:
    if f1_lines and f2_lines:
        assert len(f1_lines) == 4
        assert len(f2_lines) == 4
        f_out.write(f1_lines[0])
        f_out.write(f1_lines[1])
        f_out.write(f1_lines[2])
        f_out.write(f1_lines[3])
        f_out.write(f2_lines[0])
        f_out.write(f2_lines[1])
        f_out.write(f2_lines[2])
        f_out.write(f2_lines[3])

        f1_lines = []
        f2_lines = []

    return f1_lines, f2_lines


def join_interleaved(
    f_out: TextIO,
    f1: TextIO,
    f2: TextIO,
    strict: bool = False,
) -> None:
    readlines = True
    ix = 0
    f1_lines: list[str] = []
    f2_lines: list[str] = []

    while readlines:
        f1_line = f1.readline()
        f2_line = f2.readline()

        if f1_line == "":
            readlines = False
            if f2_line != "":
                raise Exception("Input FASTQ files do not match in length.")
            break

        if ix % 4 == 0:  # Header #1
            if strict:
                enforce_headers(f1_line, f2_line)

            f1_lines, f2_lines = flush_buffer(f_out, f1_lines, f2_lines)

        # Fill buffer up to 4 lines
        f1_lines.append(f1_line)
        f2_lines.append(f2_line)
        ix += 1

    flush_buffer(f_out, f1_lines, f2_lines)

    f1.close()
    f2.close()
    f_out.close()


def open_fastq(path: str, mode: str, use_gzip: bool = False) -> TextIO:
    if use_gzip or ".gz" in path or ".gzip" in path:
        return cast(TextIO, gzip.open(path, mode=mode))
    return open(path, mode=mode)


@click.command()
@click.argument("fastq1", type=click.Path(exists=True))
@click.argument("fastq2", type=click.Path(exists=True))
@click.argument("output_fastq", type=click.Path(), required=False, default="/dev/stdout")
@click.option(
    "--strict",
    is_flag=True,
    help="Enforce that the headers of the input FASTQ files match until the first whitespace delimiter (e.g., a space).",
)
@click.option("--gzip", "use_gzip", is_flag=True, help="Gzip output_fastq.")
def cli(fastq1: str, fastq2: str, output_fastq: str, strict: bool, use_gzip: bool):
    """Join two FASTQ paired end read files. Defaults to interleaving reads.

    Note: Requires input files to be sorted.

    FASTQ1: First input FASTQ.
    FASTQ2: Second input FASTQ.
    OUTPUT_FASTQ: Output FASTQ file name (optional, streams to STDOUT if missing).
    """

    f1 = open_fastq(fastq1, mode="rt")
    f2 = open_fastq(fastq2, mode="rt")
    f_out = open_fastq(output_fastq, mode="wt", use_gzip=use_gzip)

    join_interleaved(f_out=f_out, f1=f1, f2=f2, strict=strict)
