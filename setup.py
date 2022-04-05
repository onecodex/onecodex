#!/usr/bin/env python
"""
``onecodex``
------------

``onecodex`` provides a command line client and Python library for
interacting with the One Codex API.



Links
`````
* `One Codex: <https://www.onecodex.com/>`
* `API Docs: <http://docs.onecodex.com/>`

"""
from setuptools import setup, find_packages

with open("onecodex/version.py") as import_file:
    exec(import_file.read())


with open("README.md") as readme:
    README = readme.read()


# Dependencies
TESTING_DEPS = [
    "coverage",
    "codecov",
    "pytest",
    "pytest-cov",
    "responses",
    "testfixtures",
    "mock",
    "vega_datasets",
    "pre-commit",
    "pdfplumber",
]
ALL_DEPS = [
    "altair==4.2.0",
    "numpy>=1.11.0",
    "pandas>=1.0.3",
    "pillow>=9.0.1",
    "scikit-bio>=0.5.7",
    "scikit-learn>=0.19.0",
    "jupyter-client==7.2.0",
]
REPORT_DEPS = ALL_DEPS + [
    "notebook==6.4.10",
    "nbconvert>=6",
    # Hardcoding jinja2 version for https://github.com/jupyter/nbconvert/issues/1742
    "jinja2<=3.1",
    "WeasyPrint==54.2",
    "altair_saver==0.5.0",
]

setup(
    name="onecodex",
    version=__version__,  # noqa
    packages=find_packages(exclude=["*test*"]),
    python_requires=">=3.7",
    install_requires=[
        "boto3>=1.17.98",
        "click>=8.0",
        "jsonschema>=2.4",
        "python-dateutil>=2.5.3",
        "pytz>=2014.1",
        "sentry-sdk>=0.10.2",
        "requests>=2.27.1",
        "requests_toolbelt>=0.7.0",
        "six>=1.10.0",
        "unidecode>=1.0.23",
        "filelock>=3.0.12,<4",
    ],
    include_package_data=True,
    zip_safe=False,
    extras_require={
        "all": ALL_DEPS,
        "reports": REPORT_DEPS,
        "testing": TESTING_DEPS,
    },
    author="One Codex",
    author_email="opensource@onecodex.com",
    description="One Codex API client and Python library",
    long_description=README,
    long_description_content_type="text/markdown",
    license="MIT License",
    keywords="One Codex API Client",
    url="https://github.com/onecodex/onecodex",
    classifiers=[
        "Environment :: Console",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Internet :: WWW/HTTP",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    entry_points={
        "console_scripts": ["onecodex = onecodex.cli:onecodex"],
        "nbconvert.exporters": [
            "onecodex_pdf = onecodex.notebooks.exporters:OneCodexPDFExporter",
            "onecodex_html = onecodex.notebooks.exporters:OneCodexHTMLExporter",
            "onecodex_doc = onecodex.notebooks.exporters:OneCodexDocumentExporter",
        ],
    },
    test_suite="tests",
)
