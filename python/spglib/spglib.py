"""Python bindings for C library for finding and handling crystal."""
# Copyright (C) 2015 Atsushi Togo
# This file is part of spglib.
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations

import dataclasses
import os
import warnings
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Union

import numpy as np

import spglib

from . import _spglib
from ._compat.typing import TypeAlias
from ._compat.warnings import deprecated
from .error import SpglibError

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Any

    from numpy.typing import ArrayLike

__all__ = [
    "find_primitive",
    "get_error_message",
    "get_hall_number_from_symmetry",
    "get_version",
    "refine_cell",
    "spg_get_commit",
    "spg_get_version",
    "spg_get_version_full",
    "standardize_cell",
]

warnings.filterwarnings(
    "module", category=DeprecationWarning, message=r"dict interface.*"
)

Lattice: TypeAlias = Sequence[Sequence[float]]
Positions: TypeAlias = Sequence[Sequence[float]]
Numbers: TypeAlias = Sequence[int]
Magmoms: TypeAlias = Union[Sequence[float], Sequence[Sequence[float]]]
Cell: TypeAlias = Union["spglib.spg.SpgCell", "spglib.msg.MsgCell"]
"""Either SpgCell or MsgCell."""


OLD_ERROR_HANDLING: bool = True
"""
Use the old error handling.

Note that this variable may be removed in the future or change value in the future.
You can also use :envvar:`SPGLIB_OLD_ERROR_HANDLING` instead of altering this value.
"""


_spglib_error = ""


@dataclasses.dataclass(eq=True, frozen=True)
class DictInterface(Mapping[str, "Any"]):
    """Base class for dataclass with dict interface.

    .. versionadded:: 2.5.0
    .. deprecated:: 2.5.0
        Dict-like interface (``obj['field']``) are deprecated.
        Please use attribute interface instead (``obj.field``)
    """

    @deprecated("dict interface is deprecated. Use attribute interface instead")
    def __getitem__(self, key: str) -> Any:
        """Return the value of the key."""
        return dataclasses.asdict(self)[key]

    def __len__(self) -> int:
        """Return the number of fields."""
        return len(dataclasses.fields(self))

    def __iter__(self) -> Iterator[str]:
        """Return an iterator over the keys."""
        return iter(dataclasses.asdict(self))


@deprecated("Use __version__ or spg_get_version instead")
def get_version() -> tuple[int, int, int]:
    """Return version number of spglib with tuple of three numbers.

    .. versionadded:: 1.8.3
    .. deprecated:: 2.3.0
        Use :py:func:`spg_get_version` and ``spglib.__version__`` instead
    """
    _set_no_error()
    return _spglib.version_tuple()


def spg_get_version() -> str:
    """Get the X.Y.Z version of the detected spglib C library.

    .. versionadded:: 2.3.0

    :return: version string
    """
    _set_no_error()
    return _spglib.version_string()


def spg_get_version_full() -> str:
    """Get the full version of the detected spglib C library.

    .. versionadded:: 2.3.0

    :return: full version string
    """
    _set_no_error()
    return _spglib.version_full()


def spg_get_commit() -> str:
    """Get the commit of the detected spglib C library.

    .. versionadded:: 2.3.0

    :return: commit string
    """
    _set_no_error()
    return _spglib.commit()


def standardize_cell(
    cell: spglib.spg.SpgCell,
    to_primitive: bool = False,
    no_idealize: bool = False,
    symprec: float = 1e-5,
    angle_tolerance: float = -1.0,
) -> spglib.spg.SpgCell | None:
    """Return standardized cell. When the search failed, ``None`` is returned.

    Parameters
    ----------
    cell, symprec, angle_tolerance:
        See the docstring of get_symmetry.
    to_primitive : bool
        If True, the standardized primitive cell is created.
    no_idealize : bool
        If True, it is disabled to idealize lengths and angles of basis vectors
        and positions of atoms according to crystal symmetry.

    Returns
    -------
    The standardized unit cell or primitive cell is returned by a tuple of
    (lattice, positions, numbers). If it fails, None is returned.

    Notes
    -----
    .. versionadded:: 1.8

    Now :func:`refine_cell` and :func:`find_primitive` are shorthands of
    this method with combinations of these options.
    About the default choice of the setting, see the documentation of ``hall_number``
    argument of :func:`get_symmetry_dataset`. More detailed explanation is
    shown in the spglib (C-API) document.

    """
    _set_no_error()

    lattice, _positions, _numbers, _ = _expand_cell(cell)

    # Atomic positions have to be specified by scaled positions for spglib.
    num_atom = len(_positions)
    positions = np.zeros((num_atom * 4, 3), dtype="double", order="C")
    positions[:num_atom] = _positions
    numbers = np.zeros(num_atom * 4, dtype="intc")
    numbers[:num_atom] = _numbers
    try:
        num_atom_std = _spglib.standardize_cell(
            lattice,
            positions,
            numbers,
            num_atom,
            int(to_primitive * 1),
            int(no_idealize * 1),
            float(symprec),
            float(angle_tolerance),
        )
    except Exception as exc:
        _set_or_throw_error(exc)
        return None
    return (
        np.array(lattice.T, dtype="double", order="C"),
        np.array(positions[:num_atom_std], dtype="double", order="C"),
        np.array(numbers[:num_atom_std], dtype="intc"),
    )


def refine_cell(
    cell: spglib.spg.SpgCell, symprec: float = 1e-5, angle_tolerance: float = -1.0
) -> spglib.spg.SpgCell | None:
    """Return refined cell. When the search failed, ``None`` is returned.

    The standardized unit cell is returned by a tuple of
    (lattice, positions, numbers).

    Notes
    -----
    .. versionchanged:: 1.8

    The detailed control of standardization of unit cell can be done using
    :func:`standardize_cell`.

    """
    _set_no_error()

    lattice, _positions, _numbers, _ = _expand_cell(cell)

    # Atomic positions have to be specified by scaled positions for spglib.
    num_atom = len(_positions)
    positions = np.zeros((num_atom * 4, 3), dtype="double", order="C")
    positions[:num_atom] = _positions
    numbers = np.zeros(num_atom * 4, dtype="intc")
    numbers[:num_atom] = _numbers
    try:
        num_atom_std = _spglib.refine_cell(
            lattice,
            positions,
            numbers,
            num_atom,
            float(symprec),
            float(angle_tolerance),
        )
    except Exception as exc:
        _set_or_throw_error(exc)
        return None
    return (
        np.array(lattice.T, dtype="double", order="C"),
        np.array(positions[:num_atom_std], dtype="double", order="C"),
        np.array(numbers[:num_atom_std], dtype="intc"),
    )


def find_primitive(
    cell: spglib.spg.SpgCell, symprec: float = 1e-5, angle_tolerance: float = -1.0
) -> spglib.spg.SpgCell | None:
    """Primitive cell is searched in the input cell. If it fails, ``None`` is returned.

    The primitive cell is returned by a tuple of (lattice, positions, numbers).

    Notes
    -----
    .. versionchanged:: 1.8

    The detailed control of standardization of unit cell can be done using
    :func:`standardize_cell`.

    """
    _set_no_error()

    lattice, positions, numbers, _ = _expand_cell(cell)

    try:
        num_atom_prim = _spglib.primitive(
            lattice, positions, numbers, float(symprec), float(angle_tolerance)
        )
    except Exception as exc:
        _set_or_throw_error(exc)
        return None
    return (
        np.array(lattice.T, dtype="double", order="C"),
        np.array(positions[:num_atom_prim], dtype="double", order="C"),
        np.array(numbers[:num_atom_prim], dtype="intc"),
    )


@deprecated("Set OLD_ERROR_HANDLING to false and catch the errors directly")
def get_error_message() -> str:
    """Return error message why spglib failed.

    .. warning::
        This method is not thread safe, i.e., only safely usable
        when calling one spglib method per process.

    Notes
    -----
    .. versionadded:: 1.9.5
    .. deprecated:: 2.7.0

    """
    return _spglib_error


def _expand_cell(
    cell: Cell,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    try:
        lattice = np.array(np.transpose(cell[0]), dtype="double", order="C")
        positions = np.array(cell[1], dtype="double", order="C")
        numbers = np.array(cell[2], dtype="intc")
        if len(cell) == 4:
            magmoms = np.array(cell[3], order="C", dtype="double")
        elif len(cell) == 3:
            magmoms = None
        else:
            raise TypeError("cell has to be a tuple of 3 or 4 elements.")

        # Sanity check
        if lattice.shape != (3, 3):
            raise TypeError("lattice has to be a (3, 3) array.")
        if not (positions.ndim == 2 and positions.shape[1] == 3):
            raise TypeError("positions has to be a (num_atoms, 3) array.")
        num_atoms = positions.shape[0]
        if numbers.ndim != 1:
            raise TypeError("numbers has to be a (num_atoms,) array.")
        if len(numbers) != num_atoms:
            raise TypeError(
                "numbers has to have the same number of atoms as positions."
            )
        if magmoms is not None:
            if len(magmoms) != num_atoms:
                raise TypeError(
                    "magmoms has to have the same number of atoms as positions."
                )
            if magmoms.ndim == 1:
                # collinear
                pass
            elif magmoms.ndim == 2:
                # non-collinear
                if magmoms.shape[1] != 3:
                    raise TypeError(
                        "non-collinear magmoms has to be a (num_atoms, 3) array."
                    )
            else:
                raise TypeError("magmoms has to be a 1D or 2D array.")
    except Exception as exc:
        # Note: these will eventually be moved to the C++ side
        # For now we just recast them to SpglibError
        raise SpglibError(f"Generic Spglib error:\n{exc}") from exc

    return (lattice, positions, numbers, magmoms)


def _check_OLD_ERROR_HANDLING() -> bool:
    env_var = os.environ.get("SPGLIB_OLD_ERROR_HANLDING")
    if env_var is not None:
        if env_var.lower() in ("false", "0"):
            return False
        return True
    return OLD_ERROR_HANDLING


def _set_or_throw_error(exc: Exception, _throw: bool = False) -> None:
    if _throw or not _check_OLD_ERROR_HANDLING():
        if isinstance(exc, SpglibError):
            # Our native errors we pass transparently
            raise exc
        # Otherwise we try to recast them to SplibError
        raise SpglibError(f"Generic Spglib error:\n{exc}") from exc
    warnings.warn(
        "Set OLD_ERROR_HANDLING to false and catch the errors directly.",
        DeprecationWarning,
        stacklevel=2,
    )
    _spglib_error = str(exc)


def _set_no_error(_throw: bool = False) -> None:
    if _throw or not _check_OLD_ERROR_HANDLING():
        return
    warnings.warn(
        "Set OLD_ERROR_HANDLING to false and catch the errors directly.",
        DeprecationWarning,
        stacklevel=2,
    )
    _spglib_error = "no error"


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
