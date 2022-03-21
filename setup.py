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
    "flake8",
    "pydocstyle",
    "pytest",
    "pytest-cov",
    "responses",
    "testfixtures",
    "mock",
    "vega_datasets",
    "black",
    "pre-commit",
    "PyPDF2",
]
ALL_DEPS = [
    "altair==4.1.0",
    "numpy>=1.11.0",
    "pandas>=1.0.3",
    "pillow>=9.0.1",
    "scikit-bio>=0.5.6",
    "scikit-learn>=0.19.0",
    "jupyter-client==6.1.12",
]
REPORT_DEPS = ALL_DEPS + ["notebook==6.0.3", "nbconvert<6", "WeasyPrint==51", "altair_saver==0.5.0"]


setup(
    name="onecodex",
    version=__version__,  # noqa
    packages=find_packages(exclude=["*test*"]),
    install_requires=[
        "boto3>=1.17.98",
        "click>=7.0",
        "jsonschema>=2.4",
        "python-dateutil>=2.5.3",
        "pytz>=2014.1",
        "sentry-sdk>=0.10.2",
        "requests>=2.27.1",
        "requests_toolbelt>=0.7.0",
        "six>=1.10.0",
        "unidecode==1.0.23",
        "filelock>=3.0.12,<4",
    ],
    include_package_data=True,
    zip_safe=False,
    extras_require={
        ':python_version == "2.7"': ["futures", "enum34"],
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
