"""setup.py"""
# -*- coding: utf-8 -*-
import os
import re

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

version_file = os.path.join('moca', 'version.py')
fversion = None
metadata = None

with open(version_file, 'r') as f:
    metadata = dict(re.findall("__([a-z]+)__ = '([^']+)'", f.read()))

assert metadata is not None
fversion = metadata['version'].split('.')

test_requirements = ['pytest', 'pytest-cov', 'pytest-mpl']

MAJOR      = fversion[0]
MINOR      = fversion[1]
MICRO      = fversion[2]
ISRELEASED = False
VERSION    = '%s.%s.%s' % (MAJOR, MINOR, MICRO)

setup(
    name='moca',
    version=VERSION,
    description="Tool for motif conservation analysis",
    long_description=readme + '\n\n' + history,
    author="Saket Choudhary",
    author_email='saketkc@gmail.com',
    url='https://github.com/saketkc/moca',
    packages=[
        'moca',
        'moca.helpers',
        'moca.bedoperations',
        'moca.wigoperations',
        'moca.pipeline',
        'moca.plotter',
        'scripts'
    ],
    package_dir={'moca':
                 'moca'},
    include_package_data=True,
    package_data={
        'application' : ['application.cfg.example', 'environment.yml']
    },
    install_requires=requirements,
    license="ISC",
    zip_safe=False,
    keywords='moca',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS'
    ],
    test_suite='nose2.collector.collector',
    tests_require=test_requirements,
    entry_points = '''
            [console_scripts]
            moca=scripts.mocacli:cli
    ''',
)
