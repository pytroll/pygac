The Pygac API
=============

.. inheritance-diagram:: pygac.gac_pod.GACPODReader pygac.gac_klm.GACKLMReader
                         pygac.lac_pod.LACPODReader pygac.lac_klm.LACKLMReader


Base Classes
------------

Common functionality shared by multiple readers.

Reader
~~~~~~

.. automodule:: pygac.reader
   :members:
   :undoc-members:

GAC format reader
~~~~~~~~~~~~~~~~~

.. automodule:: pygac.gac_reader
   :members:
   :undoc-members:


LAC format reader
~~~~~~~~~~~~~~~~~

.. automodule:: pygac.lac_reader
   :members:
   :undoc-members:


POD series reader
~~~~~~~~~~~~~~~~~

.. automodule:: pygac.pod_reader
   :members:
   :undoc-members:
   :exclude-members: tsm_affected_intervals


KLM series reader
~~~~~~~~~~~~~~~~~

.. automodule:: pygac.klm_reader
   :members:
   :undoc-members:
   :exclude-members: tsm_affected_intervals


Actual Reader Implementations
-----------------------------

Actual reader implementations building upon the base classes.

GAC POD reader
~~~~~~~~~~~~~~

.. automodule:: pygac.gac_pod
   :members:
   :undoc-members:


GAC KLM reader
~~~~~~~~~~~~~~

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:


LAC POD reader
~~~~~~~~~~~~~~

.. automodule:: pygac.lac_pod
   :members:
   :undoc-members:


LAC KLM reader
~~~~~~~~~~~~~~

.. automodule:: pygac.lac_klm
   :members:
   :undoc-members:
