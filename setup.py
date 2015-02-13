"""
``onecodex``
------------

``onecodex`` provides a command line client for interaction with the
One Codex API.



Links
`````
* `One Codex: <https://www.onecodex.com/>`
* `API Docs: <http://docs.onecodex.com/>`

"""
from setuptools import setup


setup(
    name='onecodex',
    version='0.0.7',
    url='https://www.onecodex.com/',
    license='MIT',
    author='Reference Genomics, Inc.',
    author_email='help@onecodex.com',
    description='One Codex Command Line Client',
    long_description=__doc__,
    packages=['onecodex'],
    zip_safe=True,
    platforms='any',
    install_requires=[
        'requests>=2.3.0',
        'requests-toolbelt>=0.3.1'
    ],
    test_suite='nose.collector',
    entry_points={
        'console_scripts': ['onecodex = onecodex.cli:main']
    },
)
