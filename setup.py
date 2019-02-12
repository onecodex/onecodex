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
from setuptools.command.install import install


class PostInstallCommand(install):
    def run(self):
        install.run(self)

        # try to enable bash completion, if possible
        import os

        paths_to_try = ['/etc/bash_completion.d', '/usr/local/etc/bash_completion.d']

        for path in paths_to_try:
            if os.access(path, os.W_OK):
                try:
                    with open(os.path.join(path, 'onecodex'), 'w') as f:
                        f.write('eval "$(_ONECODEX_COMPLETE=source onecodex)"')
                    print('Enabled bash auto-completion for onecodex')
                    return
                except Exception:
                    print('Unable to enable bash auto-completion for onecodex')


with open('onecodex/version.py') as import_file:
    exec(import_file.read())


with open('README.md') as readme:
    README = readme.read()


setup(
    name='onecodex',
    version=__version__,  # noqa
    packages=find_packages(exclude=['*test*']),
    install_requires=[
        'boto3>=1.4.2',
        'click>=6.6',
        'jsonschema>=2.4',
        'python-dateutil>=2.5.3',
        'pytz>=2014.1',
        'raven>=6.1.0',
        'requests>=2.9',
        'requests_toolbelt>=0.7.0',
        'six>=1.10.0',
        'unidecode==1.0.23',
    ],
    include_package_data=True,
    zip_safe=False,
    extras_require={
        'all': [
            'altair==2.3.0',
            'networkx>=1.11,<2.0',
            'numpy>=1.11.0',
            'pandas>=0.23.0',
            'scikit-bio>=0.4.2,<0.5.0',
            'scikit-learn>=0.19.0',
        ],
        'testing': [
            'coverage~=4.5',
            'coveralls~=1.5',
            'flake8',
            'pytest~=4.1',
            'pytest-cov~=2.6',
            'responses',
            'testfixtures',
            'tox-pyenv==1.0.3',
            'tox>=3.5.3',
            'mock==2.0.0',
        ],
    },
    dependency_links=[],
    author='Kyle McChesney & Nick Greenfield & Roderick Bovee',
    author_email='opensource@onecodex.com',
    cmdclass={'install': PostInstallCommand},
    description='One Codex API client and Python library',
    long_description=README,
    long_description_content_type='text/markdown',
    license='MIT License',
    keywords='One Codex API Client',
    url='https://github.com/onecodex/onecodex',
    classifiers=[
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Internet :: WWW/HTTP',
    ],
    entry_points={
        'console_scripts': ['onecodex = onecodex.cli:onecodex'],
        'nbconvert.exporters': [
            'onecodex_pdf = onecodex.notebooks.exporters:OneCodexPDFExporter',
            'onecodex_html = onecodex.notebooks.exporters:OneCodexHTMLExporter',
            'onecodex_doc = onecodex.notebooks.exporters:OneCodexDocumentExporter',
            ],
    },
    test_suite='tests'
)
