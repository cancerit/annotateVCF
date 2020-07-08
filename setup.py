#!/usr/bin/env python3

from setuptools import setup

config = {
    'version': '1.0.0',
    'name': 'annotateVcf',
    'description': 'Tool to annotate vcf files...',
    'author': 'Shriram Bhosle',
    'url': 'https://github.com/cancerIT/annotateVcf',
    'author_email': 'cgphelp@sanger.ac.uk',
    'python_requires': '>= 3.3',
    'setup_requires': ['pytest','pytest-cover', 'radon'],
    'install_requires': ['tzlocal'],
    'packages': ['annotate'],
    'package_data': {'annotate':['config/*.conf','config/drvData/*']},
    'entry_points': {
        'console_scripts': ['annotateVcf=annotate.annotate_cmd:main'],
    }
}

setup(**config)

