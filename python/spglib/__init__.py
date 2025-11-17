"""Python bindings for C library for finding and handling crystal."""
# Copyright (C) 2015 Atsushi Togo
# This file is part of spglib.
# SPDX-License-Identifier: BSD-3-Clause

# TODO: _internal should not be exposed, we only continue to do so for
#  backwards compatibility
from ._internal import *  # noqa: F403
from ._version import __version__, __version_tuple__  # noqa: F401
from .kpoints import *  # noqa: F403

# fmt: off
from .spglib import (  # noqa: F401
    MagneticSpaceGroupType,
    SpaceGroupType,
    SpglibDataset,
    SpglibMagneticDataset,
    delaunay_reduce,
    find_primitive,
    get_error_message,
    get_hall_number_from_symmetry,
    get_magnetic_spacegroup_type,
    get_magnetic_spacegroup_type_from_symmetry,
    get_magnetic_symmetry,
    get_magnetic_symmetry_dataset,
    get_magnetic_symmetry_from_database,
    get_spacegroup,
    get_spacegroup_type,
    get_spacegroup_type_from_symmetry,
    get_symmetry,
    get_symmetry_dataset,
    get_symmetry_from_database,
    get_version,
    niggli_reduce,
    refine_cell,
    spg_get_commit,
    spg_get_version,
    spg_get_version_full,
    standardize_cell,
)
# fmt: on
