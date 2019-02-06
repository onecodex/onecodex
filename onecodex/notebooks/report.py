import datetime

from onecodex.exceptions import OneCodexException


class set_style(object):
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


class set_center_header(object):
    def __init__(self, text, style=None):
        self.text = text
        self.style = '' if style is None else style

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<div id="centerheader" style="{}">{}</div>'.format(self.style, self.text)}


class set_date(object):
    def __init__(self, date=None, style=None):
        self.date = '{dt:%B} {dt.day}, {dt:%Y}'.format(dt=datetime.datetime.now()) if date is None else date
        self.style = '' if style is None else style

        try:
            ipy = get_ipython()
            ipy.meta['customdate'] = self.date
        except NameError:
            pass

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


class set_logo(object):
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
        return {'text/html': '<img src="{}" width="120px" class="logo {}" style="{}" />'.format(self.url, self.classes, self.style)}


class legend(object):
    def __init__(self, text, heading=None, fignum=None, style=None):
        self.heading = '' if heading is None else '{} '.format(heading)
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
            'text/html': '<div class="figurelegend" style="{}"><b>Figure {}. {}</b>{}</div>'.format(
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
        self.style = '' if style is None else style

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
        cites = [
            '<div>'
            '<h4>References</h4>',
            '<dl class="bibliography" style="{}">'.format(self.style)
        ]

        for ref_label, (ref_num, ref_text) in self.ref_list.items():
            cites.append('<dt>{}</dt>'.format(ref_num))
            cites.append('<dd>{}</dd>'.format(ref_text))

        cites.append('</dl>')
        cites.append('</div>')

        return {'text/html': '\n'.join(cites)}


class page_break(object):
    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {'text/html': '<div class="pagebreak"></div>'}


class cover_sheet(object):
    def __init__(self, title, prepared_for, prepared_by, project_details, add_br=True):
        if add_br:
            title = title.replace('\n', '<br />')
            prepared_for = prepared_for.replace('\n', '<br />')
            prepared_by = prepared_by.replace('\n', '<br />')

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
                project_details = project_details.replace('\n', '<br />')

            self.project_details = project_details

        try:
            ipy = get_ipython()
            proj_date = ipy.meta.get('customdate')
        except NameError:
            proj_date = None

        if proj_date is None:
            proj_date = set_date().date

        self.proj_date = proj_date

    def display(self):
        from IPython.display import display
        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        body = """
        <div class="pagebreak">
          <h2 class="coverpage-title">{}</h2>
          <div class="coverpage-date">{}</div>
        """

        if self.prepared_for:
            body += """
            <h4>PREPARED FOR</h4>
            {}
            <br />
            """

        if self.prepared_by:
            body += """
            <h4>PREPARED BY</h4>
            {}
            <br />
            """

        if self.project_details:
            body += """
            <h4>PROJECT DETAILS</h4>
            {}
            """

        body += "</div>"

        body = body.format(
            self.title,
            self.proj_date,
            self.prepared_for,
            self.prepared_by,
            self.project_details,
        )

        return {'text/html': body}
