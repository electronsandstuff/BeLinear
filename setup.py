#!/usr/bin/env python

'''
BeLinear - Performant numerical solutions of the paraxial ray equation
Copyright (C) 2020 Christopher M. Pierce (contact@chris-pierce.com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import setuptools

# Load my long description
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
        name='belinear',
        version='1.0',
        description='Performant numerical solutions of the paraxial ray equation',
        author='Christopher M. Pierce',
        author_email='contact@chris-pierce.com',
        python_requires='>=3.1',
        packages=setuptools.find_packages(),
        long_description=long_description,
        long_description_content_type="text/markdown",
        install_requires = [
            'numpy',
            'matprod',
            'scipy'
        ],
        package_data={'': ['tests/test_data.pickle']},
        include_package_data=True,
        license = 'GNU Affero General Public License v3 or later (AGPLv3+)',
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
            "Development Status :: 4 - Beta",
            "Operating System :: OS Independent",
        ],
    )
