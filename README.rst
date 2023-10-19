Biolearn
========

Biolearn enables easy and versatile analyses of biomarkers of aging data. It provides tools to easily load data from publicly available sources like the 
`Gene Expression Omnibus <https://www.ncbi.nlm.nih.gov/geo/>`_, `National Health and Nutrition Examimation Survey <https://www.cdc.gov/nchs/nhanes/index.htm>`_,
and the `Framingham Heart Study <https://www.framinghamheartstudy.org/>`_. Biolearn also contains reference implemenations for common aging clock such at the 
Horvath clock, DunedinPACE and many others that can easily be run in only a few lines of code.


.. warning::

    This is a prerelease version of the biolearn library. There may be bugs and interfaces are subject to change.


Important links
===============

- Source code: https://github.com/bio-learn/biolearn/
- Documentation Homepage: https://bio-learn.github.io/

Requirements
============

Python 3.10+

Install
=======
Install biolearn using pip.

.. code-block:: bash

    pip install biolearn

To verify the library was installed correctly open python or a jupyter notebook and run:

.. code-block:: python

    from biolearn.data_library import DataLibrary

If it executes with no errors then the library is installed. To get started check out `some code examples <https://bio-learn.github.io/auto_examples/index.html>`_

Discord server
==============

The biolearn team has a `discord server <https://discord.gg/345kHcsJy>`_ to answer questions,
discuss feature requests, or have any biolearn related discussions.


Contributing
============

Detailed instructions on developer setup and how to contribute are available `in the repo <https://github.com/bio-learn/biolearn/blob/master/DEVELOPMENT.md>`_