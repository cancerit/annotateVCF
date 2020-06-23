#!/usr/bin/env python3

from setuptools import setup

config = {
    'version': '1.0.0',
    'name': 'annotateVcf',
    'description': 'Tool to annotate vcf files...',
    'author': 'Shriram Bhosle',
    'url': 'https://gitlab.com/translation/annotatevcf',
    'author_email': 'cgphelp@sanger.ac.uk',
    'python_requires': '>= 3.3',
    'setup_requires': ['pytest','pytest-cover', 'radon'],
    'install_requires': ['tzlocal'],
    'packages': ['annotateVcf'],
    'package_data': {'annotateVcf':['config/*.conf','config/*.tbi','config/*.tab.gz']},
    'entry_points': {
        'console_scripts': ['annotateVcf=annotate.annotate_cmd:main'],
    }
}

setup(**config)

