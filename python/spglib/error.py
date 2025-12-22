"""Exception types raised in the python API.

.. warning::
   Exception support is still a work in progress. Not all exceptions
   may be raised as :py:exc:`SpglibError`. Please report any exceptions
   that are not caught by :py:exc:`SpglibError`.

.. deprecated:: 2.7.0
   Currently the exception raising is opt-in, requiring to set
   :py:attr:`spglib.OLD_ERROR_HANDLING` to ``False``, otherwise
   all exceptions are redirected to :py:func:`spglib.get_error_message`.

   Starting from version 2.8.0 the default option will be flipped to always
   throw these exceptions, and the old error handling will be removed in
   version 3.0. Please test out the functionality and provide feedback on the
   GitHub issues.

.. versionadded:: 2.7.0
"""

# Copyright (C) 2015 Atsushi Togo
# This file is part of spglib.
# SPDX-License-Identifier: BSD-3-Clause
from __future__ import annotations

import os
import warnings

from ._compat.warnings import deprecated

__all__ = [
    "SpglibError",
    "get_error_message",
]


class SpglibError(Exception):
    """Base exception type for all errors raised by spglib."""


# TODO: Provide more exception types.


# Old deprecated error functions
OLD_ERROR_HANDLING: bool = True
"""
Use the old error handling.

Note that this variable may be removed in the future or change value in the future.
You can also use :envvar:`SPGLIB_OLD_ERROR_HANDLING` instead of altering this value.
"""


_spglib_error = ""


def _check_OLD_ERROR_HANDLING() -> bool:
    env_var = os.environ.get("SPGLIB_OLD_ERROR_HANDLING")
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
