# Python API exceptions

This feature is available starting in v2.7. To enable it, add the following
snippet to your code:

```python
try:
    spglib.error.OLD_ERROR_HANDLING = False
except AttributeError:
    pass
```

Note that `spglib.OLD_ERROR_HANDLING` may be removed or its default value may
change in future versions. You can also set the environment variable
`SPGLIB_OLD_ERROR_HANDLING` (`0` or `false`) instead of modifying this value
directly.

```{eval-rst}
.. automodule:: spglib.error
   :members:
   :no-index:
```
