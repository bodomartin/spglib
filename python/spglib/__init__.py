"""Python bindings for C library for finding and handling crystal."""
# Copyright (C) 2015 Atsushi Togo
# This file is part of spglib.
# SPDX-License-Identifier: BSD-3-Clause

# TODO: _internal should not be exposed, we only continue to do so for
#  backwards compatibility
from ._internal import *  # noqa: F403
from ._version import __version__, __version_tuple__  # noqa: F401
from .cell import *  # noqa: F403
from .kpoints import *  # noqa: F403
from .msg import *  # noqa: F403
from .reduce import *  # noqa: F403
from .spg import *  # noqa: F403

# fmt: off
from .spglib import (  # noqa: F401
    get_error_message,
    get_hall_number_from_symmetry,
    get_version,
    spg_get_commit,
    spg_get_version,
    spg_get_version_full,
)
# fmt: on
