Metadata Standard & Quick Search
================================

Sex codes
---------

============ =========== ==========
Label        Numeric     Meaning
------------ ----------- ----------
female       0           biological XX / reported female
male         1           biological XY / reported male
unknown      NaN         missing / indeterminate
============ =========== ==========

Age
---

* ``age`` is stored as **float, years** (e.g. ``42.7``).

Python example
--------------

.. code-block:: python

   from biolearn.data_library import GeoData

   # Search for datasets without loading them
   df = GeoData.search(sex="male", min_age=65)
   print(df.head())

   # Preview available datasets
   female_datasets = GeoData.search(sex="female")
   print(f"Found {len(female_datasets)} datasets with female subjects")
