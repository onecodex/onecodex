import pytest

pytest.importorskip("builtins")  # noqa

import builtins
import json
import mock
import os
import sys
from xml.etree import ElementTree as ET

from onecodex.exceptions import OneCodexException


class mock_nbconvert(object):
    class exporters(object):
        class html(object):
            class HTMLExporter(object):
                def __init__(self, *args, **kwargs):
                    pass

                def from_notebook_node(self, nb, resources=None, **kw):
                    return nb, resources

                def register_preprocessor(self, *args, **kwargs):
                    pass


class mock_traitlets(object):
    @staticmethod
    def default(param):
        def _wrapper(fn):
            return fn

        return _wrapper


class mock_ipython(object):
    meta = {}

    class display(object):
        def display(self, obj):
            return obj._repr_mimebundle_()


get_ipython = mock_ipython


class mock_weasyprint(object):
    class HTML(object):
        def __init__(self, string=None):
            self.text = json.dumps(string.cells)

        def write_pdf(self, buf, stylesheets=None):
            buf.write(builtins.bytes(self.text, encoding="UTF-8"))

    def CSS(self, path):
        assert os.path.exists(path)


class notebook(object):
    def __init__(self):
        self.cells = [
            {
                "cell_type": "code",
                "outputs": [
                    {
                        "metadata": {"jupyter-vega": "1234567890"},
                        "data": {"image/svg+xml": "1234567890"},
                    }
                ],
            },
            {
                "cell_type": "code",
                "outputs": [
                    {
                        "metadata": {"onecodex": "head.style"},
                        "data": {"text/css": "h1 { text-align: center; }"},
                    }
                ],
            },
            {
                "cell_type": "code",
                "outputs": [
                    {
                        "metadata": {"onecodex": "customdate"},
                        "data": {
                            "text/html": '<div id="reportdate" style="">January 1, 1900</div>'
                        },
                    }
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {"variables": {"quick": "slow", "lazy": "energetic"}},
                "source": "the {{quick}} brown fox jumped over the {{lazy}} dog",
            },
        ]


@pytest.fixture(scope="function")
def mock_reports():
    # the whole tree must be mocked for import mocking to work in python 2
    fakemodules = {
        "nbconvert": mock_nbconvert(),
        "nbconvert.exporters": mock_nbconvert().exporters(),
        "nbconvert.exporters.html": mock_nbconvert().exporters().html(),
        "traitlets": mock_traitlets(),
        "IPython": mock_ipython(),
        "IPython.display": mock_ipython().display(),
        "weasyprint": mock_weasyprint(),
    }

    # reset metadata in get_ipython() between tests
    mock_ipython().meta.clear()

    with mock.patch.dict(sys.modules, fakemodules):
        # for python2 mocking of get_ipython. import must be within mocked sys.modules
        from onecodex.notebooks import report

        report.get_ipython = None

        with mock.patch("onecodex.notebooks.report.get_ipython", mock_ipython):
            yield


@pytest.fixture(scope="function")
def mock_notebook():
    return notebook()


@pytest.mark.parametrize("metaval", [1, 2])
def test_fixture(mock_reports, metaval):
    # changes to get_ipython should be consistent within test but not between tests
    assert get_ipython().meta.get("test") is None
    get_ipython().meta["test"] = metaval
    assert get_ipython().meta.get("test") == metaval


def test_set_style(mock_reports):
    from onecodex.notebooks import report

    mystyle = "h1 { text-align: center; }"
    obj = report.set_style(mystyle)
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 2
    assert len(bundle[0]) == len(bundle[1]) == 1
    assert bundle[1]["onecodex"] == "head.style"
    ET.fromstring(bundle[0]["text/css"])
    obj.display()


def test_set_center_header(mock_reports):
    from onecodex.notebooks import report

    mytext = "center header"
    mystyle = "text-align: center"
    obj = report.set_center_header(mytext, mystyle)
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 1
    assert mytext in bundle["text/html"]
    assert mystyle in bundle["text/html"]
    ET.fromstring(bundle["text/html"])
    obj.display()


def test_set_date(mock_reports):
    from onecodex.notebooks import report

    mydate = "January 1, 1900"
    mystyle = "text-align: center"
    obj = report.set_date(mydate, mystyle)
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 2
    assert mydate in bundle[0]["text/html"]
    assert mystyle in bundle[0]["text/html"]
    assert bundle[1]["onecodex"] == "customdate"
    assert get_ipython().meta["customdate"] == mydate

    ET.fromstring(bundle[0]["text/html"])
    obj.display()


def test_title(mock_reports):
    from onecodex.notebooks import report

    mytext = "my title"
    mystyle = "text-align: center"
    obj = report.title(mytext, mystyle)
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 1
    assert mytext in bundle["text/html"]
    assert mystyle in bundle["text/html"]
    ET.fromstring(bundle["text/html"])
    obj.display()


def test_set_logo(mock_reports):
    from onecodex.notebooks import report

    myurl = "file:///path/to/logo.png"
    mystyle = "text-align: center"
    obj = report.set_logo(myurl, style=mystyle)
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 1
    assert myurl in bundle["text/html"]
    assert mystyle in bundle["text/html"]
    assert "logo-left" in bundle["text/html"]
    ET.fromstring(bundle["text/html"])
    obj.display()

    report.set_logo(myurl, position="left")
    report.set_logo(myurl, position="right")
    report.set_logo(myurl, position="center")

    with pytest.raises(OneCodexException) as e:
        report.set_logo(myurl, position="bottom")
    assert "position must be one of" in str(e.value)


def test_legend(mock_reports):
    from onecodex.notebooks import report

    mytext = "first legend"
    myheading = "first heading"
    mystyle = "text-align: center"
    obj = report.legend(mytext, heading=myheading, style=mystyle)
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 1
    assert mytext in bundle["text/html"]
    assert myheading in bundle["text/html"]
    assert mystyle in bundle["text/html"]
    ET.fromstring(bundle["text/html"])
    obj.display()

    assert get_ipython().meta["figure_count"] == 1
    report.legend("second legend", "second heading")
    assert get_ipython().meta["figure_count"] == 2


def test_reference_and_biblio(mock_reports):
    from onecodex.notebooks import report

    with pytest.raises(OneCodexException) as e:
        report.reference()
    assert "at least one of" in str(e.value)

    # duplicate references don't get new  numbers
    ref1 = report.reference("wikipedia reference 1")
    assert 1 == ref1.ref_num
    ref2 = report.reference("wikipedia reference 1")
    assert 1 == ref2.ref_num

    # new references get new numbers
    ref3 = report.reference("wikipedia reference 2")
    assert 2 == ref3.ref_num

    # lookup by tags works
    ref4 = report.reference("wikipedia reference 3", "wiki")
    assert 3 == ref4.ref_num
    ref5 = report.reference(label="wiki")
    assert 3 == ref5.ref_num

    with pytest.raises(OneCodexException) as e:
        report.reference(label="wiki2")
    assert "Cannot find citation" in str(e.value)

    with pytest.raises(OneCodexException) as e:
        report.reference("wikipedia reference 4", "wiki")
    assert "already in use" in str(e.value)

    with pytest.raises(OneCodexException) as e:
        report.reference("wikipedia reference 1", "wiki1")
    assert "already in use with label=1" in str(e.value)

    mystyle = "font-family: Monospace"
    obj = report.bibliography(style=mystyle)
    bundle = obj._repr_mimebundle_()
    assert len(bundle) == 1
    assert "wikipedia reference 1" in bundle["text/html"]
    assert "wikipedia reference 2" in bundle["text/html"]
    assert "wikipedia reference 3" in bundle["text/html"]
    assert "wikipedia reference 4" not in bundle["text/html"]
    assert mystyle in bundle["text/html"]
    ET.fromstring(bundle["text/html"])
    obj.display()


def test_page_break(mock_reports):
    from onecodex.notebooks import report

    obj = report.page_break()
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 1
    assert "pagebreak" in bundle["text/html"]
    ET.fromstring(bundle["text/html"])
    obj.display()


@pytest.mark.parametrize("add_br", [False, True])
def test_cover_sheet(mock_reports, add_br):
    from onecodex.notebooks import report

    mytitle = "report title"
    prepared_for = "darth@vader.com"
    prepared_by = "luke@skywalker.com"
    project_details = {"light saber": "blue", "force": "strong"}
    obj = report.cover_sheet(mytitle, prepared_for, prepared_by, project_details, add_br=add_br)
    bundle = obj._repr_mimebundle_()

    assert len(bundle) == 1
    assert mytitle in bundle["text/html"]
    assert prepared_for in bundle["text/html"]
    assert prepared_by in bundle["text/html"]

    for k, v in project_details.items():
        assert k in bundle["text/html"]
        assert v in bundle["text/html"]

    ET.fromstring(bundle["text/html"])
    obj.display()


@pytest.mark.xfail(reason="Need to update test env")
def test_html_export(mock_reports, mock_notebook):
    from onecodex.notebooks import exporters

    obj = exporters.OneCodexHTMLExporter()
    output, resources = obj.from_notebook_node(mock_notebook)

    # svg gets encoded for embedding in HTML
    assert (
        output.cells[0]["outputs"][0]["data"]["text/html"]
        == '<img src="data:image/svg+xml;charset=utf-8;base64,MTIzNDU2Nzg5MA==">'
    )

    # style() output gets deleted and moved to <head> in HTML output
    assert output.cells[1]["outputs"][0]["data"] == {"text/plain": ""}

    # this test will only work if nbconvert and jupyter_contrib_nbextensions are installed
    # and nbconvert.exporters.HTMLExporter is /actually/ used and not mocked
    # # text in markdown blocks gets replaced with variables in metadata
    # assert output.cells[3]['source'] == 'the slow brown fox jumped over the energetic dog'

    head_block = resources["metadata"]["head_block"]

    # custom date caused suppression of today's date in header
    assert 'id="reportdate"' not in head_block

    # one codex logo is on by default. this will have to be updated if logo is changed
    assert "iVBORw0KGgoAAA" in head_block

    # custom CSS injected into head block
    assert '<style type="text/css">h1 { text-align: center; }</style>' in head_block

    # check at the very least that it's parsable
    ET.fromstring("<root>" + head_block + "</root>")

    # do it one more time with today's date and suppressing one codex logo
    del mock_notebook.cells[2]
    patched_env = os.environ.copy()
    patched_env["ONE_CODEX_REPORT_NO_LOGO"] = "True"

    with mock.patch.object(os, "environ", patched_env):
        output, resources = obj.from_notebook_node(mock_notebook)
        head_block = resources["metadata"]["head_block"]

        assert "reportdate" in head_block
        assert "one_codex_logo.png" not in head_block


@pytest.mark.xfail(reason="Need to update test env")
def test_pdf_export(mock_reports, mock_notebook):
    from onecodex.notebooks import exporters

    obj = exporters.OneCodexPDFExporter()
    output, resources = obj.from_notebook_node(mock_notebook)

    # not much to do here without actually importing weasyprint
    assert len(output) == 605


def test_doc_portal_export():
    # TODO: test document portal export after class is implemented in reports.py
    pass
