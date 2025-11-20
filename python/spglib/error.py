"""Exception types raised in the python API.

.. warning::
   Exception support is still a work in progress. Not all exceptions
   may be raised as :py:exc:`SpglibError`. Please report any exceptions
   that are not caught by :py:exc:`SpglibError`.

.. deprecated:: 2.7.0
   Currently the exception raising is opt-in, requiring to set
   :py:attr:`spglib.spglib.OLD_ERROR_HANDLING` to ``False``, otherwise
   all exceptions are redirected to :py:func:`spglib.spglib.get_error_message`.

   Starting from version 2.8.0 the default option will be flipped to always
   throw these exceptions, and the old error handling will be removed in
   version 3.0. Please test out the functionality and provide feedback on the
   GitHub issues.

.. versionadded:: 2.7.0
"""


class SpglibError(Exception):
    """Base exception type for all errors raised by spglib."""


# TODO: Provide more exception types.
