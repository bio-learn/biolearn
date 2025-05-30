Metadata Standard & Quick Search
================================

Sex codes
---------

============ =========== ==========
Label        Numeric int Meaning
------------ ----------- ----------
female       0           biological XX / reported female
male         1           biological XY / reported male
unknown      –1          missing / indeterminate
============ =========== ==========

Age
---

* ``age`` is stored as **float, years** (e.g. ``42.7``).

Command-line example
--------------------

.. code-block:: bash

   biolearn search-metadata --field sex --value female --min-age 70

Python example
--------------

.. code-block:: python

   from biolearn.metadata import search_metadata
   df = search_metadata(sex="male", min_age=65)
   print(df.head())