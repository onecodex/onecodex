import io
import pytest

import mock
import os
import nbformat
import datetime

from traitlets.config import Config

import pdfplumber
from onecodex.notebooks.exporters import HTMLExporter, OneCodexPDFExporter


@pytest.fixture
def nb_config():
    config = Config()
    config.OneCodexPDFExporter.preprocessors = ["nbconvert.preprocessors.ExecutePreprocessor"]
    config.HTMLExporter.preprocessors = ["nbconvert.preprocessors.ExecutePreprocessor"]
    yield config


@pytest.fixture
def nb():
    yield nbformat.read("notebook_examples/example.ipynb", as_version=4)


@pytest.fixture
def generate_pdf_report(capsys, nb, nb_config):
    def _generate_pdf_report():
        exporter = OneCodexPDFExporter(config=nb_config)
        body, _ = exporter.from_notebook_node(nb)
        return body

    return _generate_pdf_report


def test_html_report_generation(capsys, nb, nb_config):
    exporter = HTMLExporter(config=nb_config)
    body, resource = exporter.from_notebook_node(nb)

    assert "Traceback" not in body  # make sure there were no Python exceptions

    # Check visualization is generated
    assert "vega-embed" in body
    assert '<div id="altair-viz-' in body

    # Check no stderr from vega CLI or similar
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("insert_date", [True, False])
def test_pdf_report_generation_do_not_insert_date(generate_pdf_report, insert_date):
    patched_env = os.environ.copy()
    patched_env["ONE_CODEX_INSERT_DATE"] = str(insert_date)

    with mock.patch.object(os, "environ", patched_env):
        body = generate_pdf_report()

    pdf = pdfplumber.open(io.BytesIO(body))

    page = pdf.pages[0]
    pdf_text = page.extract_text()

    assert "Traceback" not in pdf_text  # make sure there were no Python exceptions

    timestamp = datetime.date.today().strftime("%B %-d, %Y")

    if insert_date:
        assert timestamp in pdf_text
    else:
        assert timestamp not in pdf_text


def test_pdf_report_generation(generate_pdf_report, capsys):
    patched_env = os.environ.copy()
    patched_env["ONE_CODEX_INSERT_DATE"] = "False"
    with mock.patch.object(os, "environ", patched_env):
        body = generate_pdf_report()
    pdf = pdfplumber.open(io.BytesIO(body))
    page = pdf.pages[0]
    pdf_text = page.extract_text()

    # Check that we have an image (hack is to just check file size)
    assert len(body) > 20000
    assert len(page.images) > 0

    # Check text
    assert "Example Report" in pdf_text
    assert "NOT FOR DIAGNOSTIC USE" in pdf_text
    assert "Traceback" not in pdf_text  # make sure there were no Python exceptions

    # Check no stderr from vega CLI or similar
    assert capsys.readouterr().err == ""

    # copy test report to $pwd so that it can be uploaded to github as an
    # artifact
    with open("test-report-generated.pdf", "wb") as handle:
        handle.write(body)
