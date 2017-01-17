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
    install_requires=['potion-client==2.4.2', 'requests>=2.9', 'click>=6.6',
                      'requests_toolbelt==0.7.0', 'python-dateutil>=2.5.3',
                      'six>=1.10.0', 'boto3>=1.4.2'],
    include_package_data=True,
    zip_safe=False,
    extras_require={
        'all': ['numpy>=1.11.0', 'pandas>=0.18.1', 'matplotlib>1.5.1', 'networkx>=1.11']
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
        'console_scripts': ['onecodex = onecodex.cli:onecodex']
    },
    test_suite='tests'
)
