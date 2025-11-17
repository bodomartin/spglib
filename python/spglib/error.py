"""Exception types raised in the python API.

.. warning::
   Exception support is still a work in progress. Not all exceptions
   may be raised as :py:exc:`SpglibError`. Please report any exceptions
   that are not caught by :py:exc:`SpglibError`.
"""


class SpglibError(Exception):
    """Base exception type for all errors raised by spglib."""


# TODO: Provide more exception types.
