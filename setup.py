#!/usr/bin/env python3

#
# Copyright (c) 2021
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of annotatevcf.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.

from setuptools import setup

config = {
    'version': '1.2.5',
    'name': 'annotateVcf',
    'description': 'Tool to annotate and filter vcf files...',
    'author': 'Shriram Bhosle',
    'url': 'https://github.com/cancerIT/annotateVcf',
    'author_email': 'cgphelp@sanger.ac.uk',
    'python_requires': '>= 3.9',
    'install_requires': ['tzlocal'],
    'packages': ['annotate'],
    'package_data': {'annotate':['config/*.conf','config/*.header','config/*.json']},
    'entry_points': {
        'console_scripts': ['annotateVcf=annotate.commandline:main'],
    }
}

setup(**config)

