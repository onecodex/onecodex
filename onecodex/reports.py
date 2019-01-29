from base64 import b64encode
from builtins import bytes
import copy
import datetime
from io import BytesIO
from nbconvert.exporters.html import HTMLExporter
import os
from traitlets import default
from weasyprint import HTML, CSS

ASSETS_PATH = os.path.join(os.path.dirname(__file__), 'assets')
HTML_TEMPLATE_FILE = 'notebook_template.tpl'
CSS_TEMPLATE_FILE = 'notebook_template.css'


class title(object):
    def __init__(self, text):
        self.text = text

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<h2 class="title">{}</h2>'.format(self.text)}


class legend(object):
    def __init__(self, heading, text):
        self.heading = heading
        self.text = text

        try:
            ipy = get_ipython()
            self.fignum = ipy.meta.get('figure_count', 0) + 1
        except NameError:
            pass

        ipy.meta['figure_count'] = self.fignum

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {
            'text/html': '<div class="figurelegend"><b>Figure {}.&nbsp;{}</b>&nbsp;{}</div>'.format(
                self.fignum, self.heading, self.text
            )
        }


class pagebreak(object):
    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<div class="pagebreak">'}


class OneCodexHTMLExporter(HTMLExporter):
    export_from_notebook = 'One Codex HTML Report'
    template_path = [ASSETS_PATH]

    def from_notebook_node(self, nb, resources=None, **kw):
        """Uses nbconvert's HTMLExporter to generate HTML, with slight modifications.

        Notes
        -----
        This exporter will only save cells generated with Altair/Vega if they have an SVG image type
        stored with them. This data is only stored if our fork of `ipyvega` is installed--otherwise,
        they will be low-resolution PNGs, which will not be exported.
        """
        nb = copy.deepcopy(nb)

        for cell_idx, cell in enumerate(nb.cells):
            # base64 encode SVGs otherwise Weasyprint can't render them. delete other
            # types of output in jupyter-vega cells (e.g., image/png, or javascript)
            if cell['cell_type'] == 'code':
                for out_idx, out in enumerate(cell['outputs']):
                    if out.get('metadata') and out['metadata'].get('jupyter-vega'):
                        for mime_idx, mimetype in enumerate(out.get('data', [])):
                            if mimetype == 'image/svg+xml':
                                img = b64encode(bytes(out['data']['image/svg+xml'], encoding='UTF-8')).decode()
                                img = '<img src="data:image/svg+xml;charset=utf-8;base64,%s">' % (img,)
                                out['data'] = {'image/svg+xml': img}
                                break
                        else:
                            out['data'] = {}
            # replace variables in markdown blocks with their values
            elif cell['cell_type'] == 'markdown':
                md_vars = cell['metadata'].get('variables', {})

                if isinstance(md_vars, dict):
                    for k, v in md_vars.items():
                        cell['source'] = cell['source'].replace('{{' + k + '}}', v)

        # put the current date in metadata we can access from within the export template
        if resources is None:
            resources = {'metadata': {}}
        elif 'metadata' not in resources:
            resources['metadata'] = {}

        resources['metadata']['date'] = '{dt:%B} {dt.day}, {dt:%Y}'.format(dt=datetime.datetime.now())
        output, resources = super(OneCodexHTMLExporter, self).from_notebook_node(nb, resources=resources, **kw)

        return output, resources

    @default('template_file')
    def _template_file_default(self):
        return HTML_TEMPLATE_FILE

    @default('file_extension')
    def _file_extension_default(self):
        return '.html'


class PDFExporter(OneCodexHTMLExporter):
    export_from_notebook = 'One Codex PDF Report'
    template_path = [ASSETS_PATH]

    def from_notebook_node(self, nb, resources=None, **kw):
        """Takes output of OneCodexHTMLExporter and runs Weasyprint to get a PDF."""
        nb = copy.deepcopy(nb)

        output, resources = super(PDFExporter, self).from_notebook_node(nb, resources=resources, **kw)
        buf = BytesIO()
        HTML(string=output).write_pdf(
            buf, stylesheets=[CSS(os.path.join(ASSETS_PATH, CSS_TEMPLATE_FILE))]
        )
        buf.seek(0)
        return buf.read(), resources

    @default('file_extension')
    def _file_extension_default(self):
        return '.pdf'


class DocumentExporter(PDFExporter):
    export_from_notebook = 'Export to One Codex Document Portal'

    def from_notebook_node(self, nb, resources=None, **kw):
        """Takes PDF output from PDFExporter and uploads to One Codex Documents portal.
        """
        output, resources = super(DocumentExporter, self).from_notebook_node(nb, resources=resources, **kw)

        # TODO: implement storage of PDF in Documents portal

    @default('file_extension')
    def _file_extension_default(self):
        return '.pdf'
