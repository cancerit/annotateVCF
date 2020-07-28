#!/usr/bin/env python3

from setuptools import setup

config = {
    'version': '1.0.2',
    'name': 'annotateVcf',
    'description': 'Tool to annotate and filter vcf files...',
    'author': 'Shriram Bhosle',
    'url': 'https://github.com/cancerIT/annotateVcf',
    'author_email': 'cgphelp@sanger.ac.uk',
    'python_requires': '>= 3.3',
    'install_requires': ['tzlocal'],
    'packages': ['annotate'],
    'package_data': {'annotate':['config/*.conf']},
    'entry_points': {
        'console_scripts': ['annotateVcf=annotate.commandline:main'],
    }
}

setup(**config)

