Installation
------------

You can install the latest stable release of the software via the python package index (pypi)

.. code-block:: bash

  pip install pygac


TLE files
~~~~~~~~~
The pygac package requires Two-Line Element files stored per-satellite
in files with names such as TLE_noaa19.txt. The contents should be the
historical TLEs, i.e. a concatenation of just lines 1 and 2 without the
satellite name. For example

.. code-block::

    1 23455U 94089A   01122.93455091  .00000622  00000-0  36103-3 0  7210
    2 23455  99.1771 113.3063 0008405 277.6106  82.4106 14.12671703326608
    1 23455U 94089A   01326.97611660  .00000739  00000-0  42245-3 0  9806
    2 23455  99.1886 322.4670 0009980  66.2863 293.9354 14.12871991355419
    etc

These can be downloaded from CelesTrak via the `special data request form`_.

.. _special data request form:
    https://celestrak.com/NORAD/archives/request.php


Development
~~~~~~~~~~~

For development clone the repository from github and install pygac in editable mode.

.. code-block:: bash

  git clone git://github.com/pytroll/pygac
  cd pygac
  pip install -e .[dev]

It is recommended to activate pre-commit checks.

.. code-block:: bash

  pre-commit install

The test suite can be run using pytest.

.. code-block:: bash

  pytest -vs pygac/tests
