from base64 import b64encode
from builtins import bytes
import copy
import datetime
from io import BytesIO
from nbconvert.exporters.html import HTMLExporter
import os
import pytz
from traitlets import default

from onecodex.notebooks import report


ASSETS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'assets'))
HTML_TEMPLATE_FILE = 'notebook_template.tpl'
CSS_TEMPLATE_FILE = 'notebook_template.css'


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

        # setup our dictionary that's accessible from within jinja templates
        if resources is None:
            resources = {'metadata': {}}
        elif 'metadata' not in resources:
            resources['metadata'] = {}

        # iterate over cells in the notebook and transform data as necessary
        do_not_insert_date = False

        for cell in nb.cells:
            if cell['cell_type'] == 'code':
                for out in cell['outputs']:
                    # base64 encode SVGs otherwise Weasyprint can't render them. delete other
                    # types of output in jupyter-vega cells (e.g., image/png, or javascript)
                    if out.get('metadata') and out['metadata'].get('jupyter-vega'):
                        for mimetype in out.get('data', []):
                            if mimetype == 'image/svg+xml':
                                img = b64encode(bytes(out['data']['image/svg+xml'], encoding='UTF-8')).decode()
                                img = '<img src="data:image/svg+xml;charset=utf-8;base64,%s">' % (img,)
                                out['data'] = {'image/svg+xml': img}
                                break
                        else:
                            out['data'] = {}
                    # transfer text/css blocks to HTML <head> tag
                    elif out.get('metadata') and out['metadata'].get('onecodex') == 'head.style':
                        for mimetype in out.get('data', []):
                            if mimetype == 'text/css':
                                style_block = '<style type="text/css">{}</style>'.format(out['data']['text/css'])
                                head_block = resources['metadata'].get('head_block', '') + style_block
                                resources['metadata']['head_block'] = head_block
                                break

                        # we don't want this to be output as text, so clear it
                        out['data'] = {'text/plain': ''}
                    # if there's a custom date specified, don't insert it
                    elif out.get('metadata') and out['metadata'].get('onecodex') == 'customdate':
                        do_not_insert_date = True
            # replace variables in markdown blocks with their values
            elif cell['cell_type'] == 'markdown':
                md_vars = cell['metadata'].get('variables', {})

                if isinstance(md_vars, dict):
                    for k, v in md_vars.items():
                        cell['source'] = cell['source'].replace('{{' + k + '}}', v)

        # add one codex logo unless told not to
        if not os.environ.get('ONE_CODEX_REPORT_NO_LOGO', False):
            logo_path = 'file:///' + os.path.join(ASSETS_PATH, 'one_codex_logo.png').replace('\\', '/')
            logo_html = report.set_logo(logo_path, position='right')._repr_mimebundle_()['text/html']
            head_block = resources['metadata'].get('head_block', '') + logo_html
            resources['metadata']['head_block'] = head_block

        # add today's date unless told not to (i.e. a custom date was specified)
        if not do_not_insert_date:
            date_div = report.set_date()._repr_mimebundle_()[0]['text/html']
            head_block = resources['metadata'].get('head_block', '') + date_div
            resources['metadata']['head_block'] = head_block

        # add link to our custom CSS using system path
        css_path = 'file:///' + os.path.join(ASSETS_PATH, CSS_TEMPLATE_FILE).replace('\\', '/')
        css_link = '<link rel="stylesheet" href="{}" />'.format(css_path)
        head_block = resources['metadata'].get('head_block', '') + css_link
        resources['metadata']['head_block'] = head_block

        # tag this report for traceability, if run from notebook service. these will be transferred
        # to PDF metadata if the HTML output of this function is used as input for PDF generation
        meta_tags = [
            ('dcterms.created', datetime.datetime.now(pytz.utc).isoformat())
        ]

        user_uuid = os.environ.get('ONE_CODEX_USER_UUID')
        if user_uuid is not None:
            meta_tags.append(('author', 'one_codex_user_uuid_{}'.format(user_uuid)))

        nb_uuid = os.environ.get('ONE_CODEX_NOTEBOOK_UUID')
        if nb_uuid is not None:
            meta_tags.append(('author', 'one_codex_notebook_uuid_{}'.format(nb_uuid)))

        meta_html = ''

        for meta_name, meta_val in meta_tags:
            meta_html += '<meta name="{}" content="{}" />\n'.format(meta_name, meta_val)

        head_block = resources['metadata'].get('head_block', '') + meta_html
        resources['metadata']['head_block'] = head_block

        output, resources = super(OneCodexHTMLExporter, self).from_notebook_node(nb, resources=resources, **kw)

        return output, resources

    @default('template_file')
    def _template_file_default(self):
        return HTML_TEMPLATE_FILE

    @default('file_extension')
    def _file_extension_default(self):
        return '.html'


class OneCodexPDFExporter(OneCodexHTMLExporter):
    export_from_notebook = 'One Codex PDF Report'
    template_path = [ASSETS_PATH]

    def from_notebook_node(self, nb, resources=None, **kw):
        """Takes output of OneCodexHTMLExporter and runs Weasyprint to get a PDF."""
        from weasyprint import HTML, CSS

        nb = copy.deepcopy(nb)

        output, resources = super(OneCodexPDFExporter, self).from_notebook_node(nb, resources=resources, **kw)
        buf = BytesIO()
        HTML(string=output).write_pdf(
            buf, stylesheets=[CSS(os.path.join(ASSETS_PATH, CSS_TEMPLATE_FILE))]
        )
        buf.seek(0)
        return buf.read(), resources

    @default('file_extension')
    def _file_extension_default(self):
        return '.pdf'


class OneCodexDocumentExporter(OneCodexPDFExporter):
    export_from_notebook = 'Export to One Codex Document Portal'

    def from_notebook_node(self, nb, resources=None, **kw):
        """Takes PDF output from PDFExporter and uploads to One Codex Documents portal.
        """
        output, resources = super(OneCodexDocumentExporter, self).from_notebook_node(nb, resources=resources, **kw)

        # TODO: implement storage of PDF in Documents portal

    @default('file_extension')
    def _file_extension_default(self):
        return '.pdf'
