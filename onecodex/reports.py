from base64 import b64encode
from builtins import bytes
import copy
import datetime
from io import BytesIO
from nbconvert.exporters.html import HTMLExporter
import os
import pytz
from traitlets import default
from weasyprint import HTML, CSS

from onecodex.exceptions import OneCodexException


ASSETS_PATH = os.path.join(os.path.dirname(__file__), 'assets')
HTML_TEMPLATE_FILE = 'notebook_template.tpl'
CSS_TEMPLATE_FILE = 'notebook_template.css'


class style(object):
    def __init__(self, style):
        if not style.startswith('\n'):
            style = '\n' + style
        if not style.endswith('\n'):
            style += '\n'
        self.style = style

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        block = '\n<style type="text/css">{}</style>\n'.format(self.style)

        return {'text/css': block}, {'onecodex': 'head.style'}


class centerheader(object):
    def __init__(self, text, style=None):
        self.text = text
        self.style = '' if style is None else style

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<div id="centerheader" style="{}">{}</div>'.format(self.style, self.text)}


class date(object):
    def __init__(self, date=None, style=None):
        self.date = '{dt:%B} {dt.day}, {dt:%Y}'.format(dt=datetime.datetime.now()) if date is None else date
        self.style = '' if style is None else style

        try:
            ipy = get_ipython()
            ipy.meta['customdate'] = self.date
        except NameError:
            raise OneCodexException('Must be run from within IPython')

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return ({'text/html': '<div id="reportdate" style="{}">{}</div>'.format(self.style, self.date)},
                {'onecodex': 'customdate'})


class title(object):
    def __init__(self, text, style=None):
        self.text = text
        self.style = '' if style is None else style

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<h2 class="title" style="{}">{}</h2>'.format(self.style, self.text)}


class logo(object):
    def __init__(self, url, position='left', style=None):
        self.url = url
        self.style = '' if style is None else style

        if position == 'left':
            self.classes = 'logo-left'
        elif position == 'center':
            self.classes = 'logo-center'
        elif position == 'right':
            self.classes = 'logo-right'
        else:
            raise OneCodexException('position must be one of: left, right, center')

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<img src="{}" width="120px" class="logo {}" style="{}">'.format(self.url, self.classes, self.style)}


class legend(object):
    def __init__(self, text, heading=None, fignum=None, style=None):
        self.heading = '' if heading is None else '{}&nbsp;'.format(heading)
        self.text = text
        self.style = '' if style is None else style

        if fignum is None:
            try:
                ipy = get_ipython()
                self.fignum = ipy.meta.get('figure_count', 0) + 1
            except NameError:
                raise OneCodexException('Must be run from within IPython')

            ipy.meta['figure_count'] = self.fignum
        else:
            self.fignum = fignum

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {
            'text/html': '<div class="figurelegend" style="{}"><b>Figure {}.&nbsp;{}</b>{}</div>'.format(
                self.style, self.fignum, self.heading, self.text
            )
        }


def reference(text=None, label=None):
    if text is None and label is None:
        raise OneCodexException('Please specify at least one of: text, label')

    try:
        ipy = get_ipython()
        ref_list = ipy.meta.get('references', {})
    except NameError:
        raise OneCodexException('Must be run from within IPython')

    def to_html(ref_num):
        return '<sup class="reference">{}</sup>'.format(ref_num)

    if text is not None:
        # has this reference already been cited?
        for ref_label, (ref_num, ref_text) in ref_list.items():
            if text == ref_text:
                if label is not None and label != ref_label:
                    raise OneCodexException('Citation already in use with label={}'.format(ref_label))
                else:
                    break
        else:
            # reference has not been cited. is the label already in use?
            if label is not None and label in ref_list.keys():
                raise OneCodexException('Citation label={} already in use'.format(label))

            # create the citation and assign next number
            if not ref_list:
                ref_num = 1
            else:
                ref_num = max([x[0] for x in ref_list.values()]) + 1

            if label is None:
                ref_label = ref_num
            else:
                ref_label = label

            ref_list[ref_label] = (ref_num, text)
            ipy.meta['references'] = ref_list

        return to_html(ref_num)
    elif label is not None:
        if label not in ref_list.keys():
            raise OneCodexException('Cannot find citation with label={}'.format(label))

        return to_html(ref_list[label][0])


class bibliography(object):
    def __init__(self, style=None):
        try:
            ipy = get_ipython()
            ref_list = ipy.meta.get('references', {})
        except NameError:
            raise OneCodexException('Must be run from within IPython')

        self.ref_list = ref_list

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        cites = ['<h4>References</h4>', '<dl class="bibliography">']

        for ref_label, (ref_num, ref_text) in self.ref_list.items():
            cites.append('<dt>{}</dt>'.format(ref_num))
            cites.append('<dd>{}</dd>'.format(ref_text))

        cites.append('</dl>')

        return {'text/html': '\n'.join(cites)}


class pagebreak(object):
    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<div class="pagebreak">'}


class coversheet(object):
    def __init__(self, title, prepared_for, prepared_by, project_details, add_br=True):
        if add_br:
            title = title.replace('\n', '<br>')
            prepared_for = prepared_for.replace('\n', '<br>')
            prepared_by = prepared_by.replace('\n', '<br>')

        self.title = title
        self.prepared_for = prepared_for
        self.prepared_by = prepared_by

        if isinstance(project_details, dict):
            deets = ['<dl class="coverpage-details">']

            for k, v in project_details.items():
                deets.append('<dt>{}:</dt>'.format(k))
                deets.append('<dd>{}</dd>'.format(v))

            deets.append('</dl>')

            self.project_details = '\n'.join(deets)
        else:
            if add_br:
                project_details = project_details.replace('\n', '<br>')

            self.project_details = project_details

        try:
            ipy = get_ipython()
            proj_date = ipy.meta.get('customdate')
        except NameError:
            proj_date = None

        if proj_date is None:
            proj_date = date().date

        self.proj_date = proj_date

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        body = """
        <h2 class="coverpage-title">{}</h2>
        <div class="coverpage-date">{}</div>
        <h4>PREPARED FOR</h4>
        {}
        <br>
        <h4>PREPARED BY</h4>
        {}
        <br>
        <h4>PROJECT DETAILS</h4>
        {}
        """.format(
            self.title,
            self.proj_date,
            self.prepared_for,
            self.prepared_by,
            self.project_details,
        )

        pagebreak_html = pagebreak()._repr_mimebundle_()['text/html']

        return {'text/html': body + pagebreak_html}


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
                                head_block = resources['metadata'].get('head_block', '') + out['data']['text/css']
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
            logo_html = logo(logo_path, position='right')._repr_mimebundle_()['text/html']
            head_block = resources['metadata'].get('head_block', '') + logo_html
            resources['metadata']['head_block'] = head_block

        # add today's date unless told not to (i.e. a custom date was specified)
        if not do_not_insert_date:
            date_div = date()._repr_mimebundle_()[0]['text/html']
            head_block = resources['metadata'].get('head_block', '') + date_div
            resources['metadata']['head_block'] = head_block

        # add link to our custom CSS using system path
        css_path = 'file:///' + os.path.join(ASSETS_PATH, CSS_TEMPLATE_FILE).replace('\\', '/')
        css_link = '<link rel="stylesheet" href="{}">'.format(css_path)
        head_block = resources['metadata'].get('head_block', '') + css_link
        resources['metadata']['head_block'] = head_block

        # tag this report for traceability, if run from notebook service. these will be transferred
        # to PDF metadata if the HTML output of this function is used as input for PDF generation
        meta_tags = [
            ('dcterms.created', datetime.datetime.now(pytz.utc).isoformat())
        ]

        user_uuid = os.environ.get('ONE_CODEX_USER_UUID')

        if user_uuid is not None:
            nb_uuid = os.environ.get('HOSTNAME', '').split('-')

            if len(nb_uuid) == 2:
                nb_uuid = nb_uuid[1]
            else:
                nb_uuid = 'not_found'

            meta_tags.append(('author', 'one_codex_user_uuid_{}'.format(user_uuid)))
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
