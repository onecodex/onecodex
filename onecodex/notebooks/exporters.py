from base64 import b64encode
from builtins import bytes
import copy
import datetime
from io import BytesIO
import json
from nbconvert.exporters.html import HTMLExporter
from nbconvert.preprocessors import Preprocessor
import os
import pytz
from traitlets import default

from onecodex.exceptions import UploadException
from onecodex.notebooks import report

ASSETS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "assets"))
HTML_TEMPLATE_FILE = "notebook_template.tpl"
CSS_TEMPLATE_FILE = "notebook_template.css"


class AltairPreprocessor(Preprocessor):
    """Update Altair-generated SVG cells for Weasyprint-compat encoding."""

    def convert_cell(self, cell):
        """Update any cells with SVGs/PNGs."""
        if cell["cell_type"] != "code":
            return cell

        for out in cell["outputs"]:
            if "data" not in out:
                continue
            if "image/svg+xml" in out["data"]:
                img = b64encode(bytes(out["data"]["image/svg+xml"], encoding="UTF-8")).decode()
                out["data"] = {
                    "text/html": '<img src="data:image/svg+xml;charset=utf-8;base64,{}">'.format(
                        img
                    )
                }
            elif "image/png" in out["data"]:
                out["data"] = {"text/html": '<img src="{}"'.format(out["data"]["image/png"])}
        return cell

    def preprocess(self, nb, resources):
        """Preprocess each notebook cell. See base.py for details."""
        # Filter out cells that meet the conditions
        nb.cells = [self.convert_cell(cell) for cell in nb.cells]

        return nb, resources


class OneCodexHTMLExporter(HTMLExporter):
    export_from_notebook = "One Codex HTML Report"

    def __init__(self, config=None, **kw):
        super(OneCodexHTMLExporter, self).__init__(config=config, **kw)

        # this preprocessor converts {{variable}} tags in markdown blocks to their values
        # see https://jupyter-contrib-nbextensions.readthedocs.io/en/latest/nbextensions/python-
        # markdown/readme.html
        try:
            from jupyter_contrib_nbextensions.nbconvert_support import (  # noqa
                PyMarkdownPreprocessor,
            )

            self.register_preprocessor(
                "jupyter_contrib_nbextensions.nbconvert_support.PyMarkdownPreprocessor",
                enabled=True,
            )
        except ImportError:
            return

    @property
    def template_paths(self):
        # https://github.com/jupyter/nbconvert/issues/1492
        return super()._template_paths() + [ASSETS_PATH]

    def from_notebook_node(self, nb, resources=None, **kw):
        """Apply nbconvert's HTMLExporter to generate HTML, with slight modifications.

        Notes
        -----
        This exporter will only save cells generated with Altair/Vega if they have an SVG image type
        stored with them. This data is only stored if our fork of `ipyvega` is installed or the onecodex
        # renderer is used--otherwise, they will be low-resolution PNGs, which will not be exported.
        """
        nb = copy.deepcopy(nb)

        # setup our dictionary that's accessible from within jinja templates
        if resources is None:
            resources = {"metadata": {}}
        elif "metadata" not in resources:
            resources["metadata"] = {}

        # iterate over cells in the notebook and transform data as necessary
        do_not_insert_date = False

        # add one codex logo unless told not to
        if not os.environ.get("ONE_CODEX_REPORT_NO_LOGO", False):
            img = b64encode(
                bytes(open(os.path.join(ASSETS_PATH, "one_codex_logo.png"), "rb").read())
            ).decode()
            img = "data:image/png;charset=utf-8;base64,%s" % (img,)
            logo_html = report.set_logo(img, position="right")._repr_mimebundle_()["text/html"]
            head_block = resources["metadata"].get("head_block", "") + logo_html
            resources["metadata"]["head_block"] = head_block

        # add today's date unless told not to (i.e. a custom date was specified)
        if not do_not_insert_date:
            date_div = report.set_date()._repr_mimebundle_()[0]["text/html"]
            head_block = resources["metadata"].get("head_block", "") + date_div
            resources["metadata"]["head_block"] = head_block

        # add the report ID and not for diagnostic use footer
        report_info_div = '<div id="reportinfo">{}</div>'.format(
            "{} - NOT FOR DIAGNOSTIC USE".format(os.environ["ONE_CODEX_REPORT_UUID"])
            if os.environ.get("ONE_CODEX_REPORT_UUID")
            else "NOT FOR DIAGNOSTIC USE"
        )
        head_block = resources["metadata"].get("head_block", "") + report_info_div
        resources["metadata"]["head_block"] = head_block

        # embed the default CSS
        css = open(os.path.join(ASSETS_PATH, CSS_TEMPLATE_FILE), "r").read()
        css = '<style type="text/css">{}</style>'.format(css)
        head_block = resources["metadata"].get("head_block", "") + css
        resources["metadata"]["head_block"] = head_block

        # tag this report for traceability, if run from notebook service. these will be transferred
        # to PDF metadata if the HTML output of this function is used as input for PDF generation
        meta_tags = [("dcterms.created", datetime.datetime.now(pytz.utc).isoformat())]

        user_uuid = os.environ.get("ONE_CODEX_USER_UUID")
        if user_uuid is not None:
            meta_tags.append(("author", "one_codex_user_uuid_{}".format(user_uuid)))

        nb_uuid = os.environ.get("ONE_CODEX_NOTEBOOK_UUID")
        if nb_uuid is not None:
            meta_tags.append(("author", "one_codex_notebook_uuid_{}".format(nb_uuid)))

        meta_html = ""

        for meta_name, meta_val in meta_tags:
            meta_html += '<meta name="{}" content="{}" />\n'.format(meta_name, meta_val)

        head_block = resources["metadata"].get("head_block", "") + meta_html
        resources["metadata"]["head_block"] = head_block

        output, resources = super(OneCodexHTMLExporter, self).from_notebook_node(
            nb, resources=resources, **kw
        )

        return output, resources

    @default("template_file")
    def _template_file_default(self):
        return HTML_TEMPLATE_FILE

    @default("file_extension")
    def _file_extension_default(self):
        return ".html"


class OneCodexPDFExporter(OneCodexHTMLExporter):
    export_from_notebook = "One Codex PDF Report"

    def __init__(self, config=None, **kw):
        super(OneCodexPDFExporter, self).__init__(config=config, **kw)
        self.register_preprocessor(AltairPreprocessor, enabled=True)

    @property
    def template_paths(self):
        # https://github.com/jupyter/nbconvert/issues/1492
        return super()._template_paths() + [ASSETS_PATH]

    def from_notebook_node(self, nb, resources=None, **kw):
        """Take output of OneCodexHTMLExporter and run Weasyprint to get a PDF."""
        from weasyprint import HTML, CSS

        nb = copy.deepcopy(nb)

        output, resources = super(OneCodexPDFExporter, self).from_notebook_node(
            nb, resources=resources, **kw
        )
        buf = BytesIO()
        HTML(string=output, base_url=ASSETS_PATH).write_pdf(
            buf, stylesheets=[CSS(os.path.join(ASSETS_PATH, CSS_TEMPLATE_FILE))]
        )
        buf.seek(0)
        return buf.read(), resources

    @default("file_extension")
    def _file_extension_default(self):
        return ".pdf"


class OneCodexDocumentExporter(OneCodexPDFExporter):
    # make this a string and it will appear in the 'Download as...' menu
    export_from_notebook = None

    def from_notebook_node(self, nb, resources=None, **kw):
        """Take PDF output from PDFExporter and upload to One Codex Documents portal."""
        output, resources = super(OneCodexDocumentExporter, self).from_notebook_node(
            nb, resources=resources, **kw
        )

        from onecodex import Api
        from onecodex.lib.upload import _upload_document_fileobj
        import sentry_sdk

        ocx = Api()

        default_filename = "Analysis Report - {dt:%B} {dt.day}, {dt:%Y}".format(
            dt=datetime.datetime.now()
        )

        file_name = resources["metadata"].get("one_codex_doc_portal_filename", default_filename)

        try:
            document_id = _upload_document_fileobj(
                BytesIO(output), file_name, ocx.Documents._resource
            )
        except UploadException as exc:
            resp = json.dumps({"status": 500, "message": str(exc)})
            return resp, resources
        except Exception as exc:
            sentry_sdk.capture_exception(exc)
            resp = json.dumps(
                {
                    "status": 500,
                    "message": "Upload failed. Please contact support@onecodex.com for assistance.",
                }
            )
            return resp, resources

        resp = json.dumps({"status": 200, "document_id": document_id})
        return resp, resources

    @default("file_extension")
    def _file_extension_default(self):
        return ".pdf"
