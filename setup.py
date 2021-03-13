#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the
#   Bartolina Project (https://github.com/exiliadadelsur/Bartolina).
# Copyright (c) 2020 Noelia Rocío Perez and Claudio Antonio Lopez Cortez
# License: MIT
#   Full Text: https://github.com/exiliadadelsur/Bartolina/blob/master/LICENSE


# =============================================================================
# DOCS
# =============================================================================

"""This file is for distribute and install Bartolina
"""


# =============================================================================
# IMPORTS
# =============================================================================

import os
import pathlib

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup

# =============================================================================
# CONSTANTS
# =============================================================================

REQUIREMENTS = [
    "numpy",
    "astropy",
    "pandas",
]

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))


# =============================================================================
# FUNCTIONS
# =============================================================================


def do_setup():
    setup(
        name="pehuen",
        description="Corrections for the redshift distortion",
        author=["Noelia Rocío Perez", "Claudio Antonio Lopez Cortez"],
        url="https://github.com/exiliadadelsur/Bartolina",
        license="MIT",
        keywords=["space redshift", "kaiser", "finger of god", "fog"],
        packages=["pehuen"],
        install_requires=REQUIREMENTS,
    )


if __name__ == "__main__":

    do_setup()
