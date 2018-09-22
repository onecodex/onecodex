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


with open('onecodex/version.py') as import_file:
    exec(import_file.read())


setup(
    name='onecodex',
    version=__version__,  # noqa
    packages=find_packages(exclude=['*test*']),
    install_requires=['potion-client==2.5.1', 'requests>=2.9', 'click>=6.6',
                      'requests_toolbelt==0.7.0', 'python-dateutil>=2.5.3',
                      'six>=1.10.0', 'boto3>=1.4.2', 'raven>=6.1.0', 'pytz>=2014.1'],
    include_package_data=True,
    zip_safe=False,
    extras_require={
        'all': ['numpy>=1.11.0', 'pandas>=0.20.0,<0.21.0', 'matplotlib>1.5.1',
                'seaborn>=0.8', 'scikit-learn>=0.19.0', 'scikit-bio==0.4.2',
                'networkx>=1.11'],
        'testing': ['flake8', 'testfixtures', 'responses', 'coverage', 'pytest==3.0.5',
                    'mock==2.0.0', 'pytest-cov==2.4.0', 'coveralls==1.1', 'tox-pyenv==1.0.3'],
    },
    dependency_links=[],
    author='Kyle McChesney & Nick Greenfield & Roderick Bovee',
    author_email='opensource@onecodex.com',
    long_description=__doc__,
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
    },
    test_suite='tests'
)
