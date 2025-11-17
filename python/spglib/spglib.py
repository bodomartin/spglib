"""Python bindings for C library for finding and handling crystal."""
# Copyright (C) 2015 Atsushi Togo
# This file is part of spglib.
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np

import spglib

from . import _spglib
from ._compat.typing import TypeAlias
from ._compat.warnings import deprecated
from .error import _set_no_error, _set_or_throw_error

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

__all__ = [
    "get_hall_number_from_symmetry",
]

Cell: TypeAlias = Union["spglib.spg.SpgCell", "spglib.msg.MsgCell"]
"""Either SpgCell or MsgCell."""


@deprecated("Use get_spacegroup_type_from_symmetry instead")
def get_hall_number_from_symmetry(
    rotations: ArrayLike[np.intc],
    translations: ArrayLike[np.double],
    symprec: float = 1e-5,
) -> int | None:
    """Hall number is obtained from a set of symmetry operations. If fails, return None.

    .. deprecated:: 2.0
        Replaced by {func}`get_spacegroup_type_from_symmetry`.

    Return one of ``hall_number`` corresponding to a space-group type of the given
    set of symmetry operations. When multiple ``hall_number`` exist for the
    space-group type, the smallest one (the first description of the space-group
    type in International Tables for Crystallography) is chosen. The definition of
    ``hall_number`` is found at :ref:`dataset_spg_get_dataset_spacegroup_type` and
    the corresponding space-group-type information is obtained through
    {func}`get_spacegroup_type`.

    This is expected to work well for the set of symmetry operations whose
    distortion is small. The aim of making this feature is to find
    space-group-type for the set of symmetry operations given by the other
    source than spglib.

    Note that the definition of ``symprec`` is
    different from usual one, but is given in the fractional
    coordinates and so it should be small like ``1e-5``.
    """
    _set_no_error()

    r = np.array(rotations, dtype="intc", order="C")
    t = np.array(translations, dtype="double", order="C")
    try:
        hall_number = _spglib.hall_number_from_symmetry(r, t, float(symprec))
    except Exception as exc:
        _set_or_throw_error(exc)
        return None
    return hall_number
