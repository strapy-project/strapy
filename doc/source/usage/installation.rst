Installation
=================================

Dependencies
------------
`strapy` is a python package, and a CPython 3 interpreter must be installed
before it can be used. For installation guidance for python see the python
documentation (on `Windows <https://docs.python.org/3/using/windows.html>`_, or
`Unix <https://docs.python.org/3/using/unix.html>`_ systems).

`strapy` uses `numpy` and `scipy` for matrix mathematics, `sympy` for symbolic
mathematics and `pyctmm` for transfer matrix modelling. If `strapy` is installed
with pip with the exception of `pyctmm` these dependencies will be installed
automatically. `pyctmm` must be installed before `strapy`, refer to the
:code:`ctmm`/`pyctmm` documentation for installation instructions.

Installation from source
------------------------
Once `pyctmm` has been installed, from a terminal open in the `strapy` root
folder: ::

    pip install .

Testing
-------
Once strapy is installed, to test from the strapy root directory run: ::

    cd test
    python -m unittest

Tests are provided for all optical components with the exception of
`strapy.components.Dump`. Components are tested against physical scenarios with
known solutions. Although the python unittest framework is used, these tests are
not true unit tests, as all tests contain multiple components. The aim is to test
that the implemented physics is correct, with the hope that bugs in the software
will also manifest themselves as incorrect physical results.