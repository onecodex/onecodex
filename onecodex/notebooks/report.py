import datetime

from onecodex.exceptions import OneCodexException


class set_style(object):
    """Inserts a <style> block in <head>. Used to override default styling.

    Parameters
    ----------
    style : `string`
        CSS to override default styling of the entire report.
    """

    def __init__(self, style):
        if not style.startswith("\n"):
            style = "\n" + style
        if not style.endswith("\n"):
            style += "\n"
        self.style = style

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        block = '\n<style type="text/css">{}</style>\n'.format(self.style)

        return {"text/css": block}, {"onecodex": "head.style"}


class set_center_header(object):
    """Inserts text in the center of the header at the top of every page of the report.

    Parameters
    ----------
    text : `string`
        Centered text to display to the top of every page.
    style : `string`, optional
        CSS to override default styling.
    """

    def __init__(self, text, style=None):
        self.text = text
        self.style = "" if style is None else style

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {
            "text/html": '<div id="centerheader" style="{}">{}</div>'.format(self.style, self.text)
        }


class set_date(object):
    """Set report-wide date, overriding default (today).

    Parameters
    ----------
    date : `string`
        Date to use in page footers and cover pages. Default uses January 1, 1900 format.
    style : `string`, optional
        CSS to override default styling.
    """

    def __init__(self, date=None, style=None):
        self.date = (
            "{dt:%B} {dt.day}, {dt:%Y}".format(dt=datetime.datetime.now()) if date is None else date
        )
        self.style = "" if style is None else style

        try:
            ipy = get_ipython()
            ipy.meta["customdate"] = self.date
        except NameError:
            pass

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return (
            {"text/html": '<div id="reportdate" style="{}">{}</div>'.format(self.style, self.date)},
            {"onecodex": "customdate"},
        )


class title(object):
    """Insert an <h2> title block. Used for either section or report titles.

    Parameters
    ----------
    text : `string`
        Text to insert into <h2> block.
    style : `string`, optional
        CSS to override default styling.
    """

    def __init__(self, text, style=None):
        self.text = text
        self.style = "" if style is None else style

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {"text/html": '<h1 class="title" style="{}">{}</h1>'.format(self.style, self.text)}


class set_logo(object):
    """Place a custom logo at the top of every page of the report.

    Parameters
    ----------
    url : `string`
        Path to the custom report logo as a URL. If a local file, use file://.
    position : {'left', 'center', 'right'}, optional
        Where to place the custom logo. Default is the top-left corner of every page.
    style : `string`, optional
        CSS to override default styling.
    """

    def __init__(self, url, position="left", style=None):
        self.url = url
        self.style = "" if style is None else style

        if position == "left":
            self.classes = "logo-left"
        elif position == "center":
            self.classes = "logo-center"
        elif position == "right":
            self.classes = "logo-right"
        else:
            raise OneCodexException("position must be one of: left, right, center")

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {
            "text/html": '<img src="{}" width="120px" class="logo {}" style="{}" />'.format(
                self.url, self.classes, self.style
            )
        }


class legend(object):
    """Add a figure legend. Call this after generating a figure.

    Parameters
    ----------
    text : `string`
        The meat of the figure legend.
    heading : `string`, optional
        Bolded text to appear before the meat of the legend.
    fignum : `string` or `integer`, optional
        The number of this figure. If not specified, will auto-increment every time this method is
        called, starting with Figure 1.
    style : `string`, optional
        CSS to override default styling.

    Notes
    -----
    Figure legend text looks something like this:

        <b>Figure {fignum}. {heading}</b> {text}

    The figure number will be auto-incremented every time legend() is called, but can be overriden
    by passing the fignum kwarg.
    """

    def __init__(self, text, heading=None, fignum=None, style=None):
        self.heading = "" if heading is None else "{} ".format(heading)
        self.text = text
        self.style = "" if style is None else style

        if fignum is None:
            try:
                ipy = get_ipython()
                self.fignum = ipy.meta.get("figure_count", 0) + 1
            except NameError:
                raise OneCodexException("Must be run from within IPython")

            ipy.meta["figure_count"] = self.fignum
        else:
            self.fignum = fignum

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {
            "text/html": '<div class="figurelegend" style="{}"><b>Figure {}. {}</b>{}</div>'.format(
                self.style, self.fignum, self.heading, self.text
            )
        }


class reference(object):
    """Add a reference to the bibliography and insert a superscript number.

    Parameters
    ----------
    text : `string`, optional
        The complete text of the reference, e.g. Roo, et al. "How to Python." Nature, 2019.
    label : `string`, optional
        A short label to describe this reference.

    Notes
    -----
    1) Every time reference() is called, the reference number is auto-incremented. That is, the first
       time you call this, a superscript 1 is inserted. The next time you call this (with a different
       reference), a superscript 2 is inserted.

    2) This function returns HTML. It is meant to be used inside a Markdown cell in your IPython
       notebook, or concatenated together with another string that's used as input to a function
       here in the `report` module.

    Examples
    --------
    You want to insert a reference at the current position, and store it using a short label so you
    can access it again without typing the entire reference text.

        >>> reference('Roo, et al. "How to Python." Nature, 2019.', 'roo1')
        '<sup class="reference">1</sup>'

    The next time you want to insert this reference, just use the short 'roo1' label.

        >>> reference(label='roo1')
        '<sup class="reference">1</sup>'

    You want to insert a list of references in a single IPython cell, each with a short label, and
    use them all later without displaying any of the superscript reference numbers now.

        _ = reference('Roo, et al. "How to Python." Nature, 2019.', 'roo1')
        _ = reference('Roo, et al. "The Tao of Roo." Random House, 2018.', 'roo2')
        _ = reference('Roo, et al. "Roo and the Art of Database Maintenance." N/A, 2017.', 'roo3')

        ~~~ And later, in a Markdown cell in your IPython notebook ~~~

        As Roo, et al. outlined in a previous work{reference(label='roo2')}, all play and no work
        makes for a happy dog. Later, the same authors applied similar reasoning to the care of
        Burmese Pythons{reference(label='roo1')}. By comparing the care of dogs and Pythons to
        SQL databases, Roo, et al. make a compelling argument for why writing docstrings can be fun
        and not just a chore{reference(label='roo3')}.

    You want to insert a reference into a figure legend, using `report.legend`.

        report.legend(
            'As you can clearly see in the above figure, the data supports my conclusion '
            'and does not support the conclusion of my peers{reference(label='similar_paper1')}. '
            'This is most likely because I am smarter and more attractive than the authors of '
            'those other publications{reference(label='ego_and_insecurity1')}.'
        )
    """

    def __init__(self, text=None, label=None):
        if text is None and label is None:
            raise OneCodexException("Please specify at least one of: text, label")

        self.text = text or ""
        self.label = label or ""

        try:
            ipy = get_ipython()
            self.ref_list = ipy.meta.get("references", {})
        except NameError:
            raise OneCodexException("Must be run from within IPython")

        if text:
            # has this reference already been cited?
            for ref_label, (ref_num, ref_text) in self.ref_list.items():
                print(ref_label, ref_num, ref_text)
                if text == ref_text:
                    if label and label != ref_label:
                        raise OneCodexException(
                            "Citation already in use with label={}".format(ref_label)
                        )
                    else:
                        self.ref_num = ref_num
                        break
            else:
                # reference has not been cited. is the label already in use?
                if label in self.ref_list.keys():
                    raise OneCodexException("Citation label={} already in use".format(label))

                # create the citation and assign next number
                if not self.ref_list:
                    self.ref_num = 1
                else:
                    self.ref_num = max([x[0] for x in self.ref_list.values()]) + 1

                if not label:
                    ref_label = self.ref_num
                else:
                    ref_label = label

                self.ref_list[ref_label] = (self.ref_num, text)
                ipy.meta["references"] = self.ref_list

        elif label:
            if label not in self.ref_list.keys():
                raise OneCodexException("Cannot find citation with label={}".format(label))

            self.ref_num = self.ref_list[label][0]

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {"text/html": '<sup class="reference">{}</sup>'.format(self.ref_num)}

    def __str__(self):
        return self._repr_mimebundle_()["text/html"]


class bibliography(object):
    """Adds a bibliography containing all references cited using `report.reference`.

    Parameters
    ----------
    style : `string`, optional
        CSS to override default styling.
    """

    def __init__(self, style=None):
        self.style = "" if style is None else style

        try:
            ipy = get_ipython()
            ref_list = ipy.meta.get("references", {})
        except NameError:
            raise OneCodexException("Must be run from within IPython")

        self.ref_list = ref_list

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        cites = [
            "<div>" "<h4>References</h4>",
            '<dl class="bibliography" style="{}">'.format(self.style),
        ]

        for ref_label, (ref_num, ref_text) in self.ref_list.items():
            cites.append("<dt>{}</dt>".format(ref_num))
            cites.append("<dd>{}</dd>".format(ref_text))

        cites.append("</dl>")
        cites.append("</div>")

        return {"text/html": "\n".join(cites)}


class page_break(object):
    """Inserts a page break."""

    def display(self):
        from IPython.display import display

        display(self)

    def _repr_mimebundle_(self, include=None, exclude=None):
        return {"text/html": '<div class="pagebreak"></div>'}


class cover_sheet(object):
    """Inserts a full-page cover sheet containing report title and various bits of data.

    Parameters
    ----------
    title : `string`, optional
        The title of the report.
    prepared_for : `string`, optional
        The name, address, phone number, e-mail, or other information about the person or
        organization for whom this report was generated.
    prepared_by : `string`, optional
        Same as above, but for the person or organization generating this report.
    project_details : `string` or `dict`, optional
        Information bits of data related to this project. If given a dictionary, will generate a
        two-column table where keys appear in the left column, bold, and values appear in the right
        column. If given a string, will dump that string into the project details section.
    add_br : `bool`, optional
        By default, will replace newlines with HTML <br> tags. If false, won't.
    """

    def __init__(self, title, prepared_for, prepared_by, project_details, add_br=True):
        if add_br:
            title = title.replace("\n", "<br />")
            prepared_for = prepared_for.replace("\n", "<br />")
            prepared_by = prepared_by.replace("\n", "<br />")

        self.title = title
        self.prepared_for = prepared_for
        self.prepared_by = prepared_by

        if isinstance(project_details, dict):
            deets = ['<dl class="coverpage-details">']

            for k, v in project_details.items():
                deets.append("<dt>{}:</dt>".format(k))
                deets.append("<dd>{}</dd>".format(v))

            deets.append("</dl>")

            self.project_details = "\n".join(deets)
        else:
            if add_br:
                project_details = project_details.replace("\n", "<br />")

            self.project_details = project_details

        try:
            ipy = get_ipython()
            proj_date = ipy.meta.get("customdate")
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
            self.title, self.proj_date, self.prepared_for, self.prepared_by, self.project_details
        )

        return {"text/html": body}
