Biolearn
========

Biolearn enables easy and versatile analyses of biomarkers of aging data. It provides tools to easily load data from publicly available sources like the
`Gene Expression Omnibus <https://www.ncbi.nlm.nih.gov/geo/>`_, `National Health and Nutrition Examimation Survey <https://www.cdc.gov/nchs/nhanes/index.htm>`_,
and the `Framingham Heart Study <https://www.framinghamheartstudy.org/>`_. Biolearn also contains reference implemenations for common aging clock such at the
Horvath clock, DunedinPACE and many others that can easily be run in only a few lines of code. You can read more about it in our `paper <https://www.biorxiv.org/content/10.1101/2023.12.02.569722v2>`_.


.. warning::

    This is a prerelease version of the biolearn library. There may be bugs and interfaces are subject to change.


Important links
===============

- Source code: https://github.com/bio-learn/biolearn/
- Documentation Homepage: https://bio-learn.github.io/

Requirements
============

- Python 3.10 or higher
- pip (Python package installer)
- venv (for virtual environments, recommended)

System Prerequisites
====================

**Ubuntu/Debian:**

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install python3-pip python3-venv

**Fedora/RHEL/CentOS:**

.. code-block:: bash

    sudo dnf install python3-pip python3-virtualenv

**macOS:**

Python from Homebrew includes pip and venv:

.. code-block:: bash

    brew install python@3.10

**Windows:**

Download Python from `python.org <https://www.python.org/downloads/>`_. The installer includes pip and venv by default.

Install
=======

We recommend using a virtual environment:

.. code-block:: bash

    # Create and activate virtual environment
    python3 -m venv biolearn-env
    source biolearn-env/bin/activate  # Linux/macOS
    # biolearn-env\Scripts\activate    # Windows

    # Install biolearn
    pip install biolearn

**Alternative: Install without virtual environment**

.. code-block:: bash

    pip install --user biolearn

Verify Installation
===================

Test the installation:

.. code-block:: python

    from biolearn.data_library import DataLibrary

If it executes without errors, the library is installed correctly. For a comprehensive installation test:

.. code-block:: bash

    python scripts/test_installation.py

To get started, check out `some code examples <https://bio-learn.github.io/auto_examples/index.html>`_

Discord server
==============

The biolearn team has a `discord server <https://discord.gg/wZH85WRTxN>`_ to answer questions,
discuss feature requests, or have any biolearn related discussions.

Issues
======

If you find any bugs with biolearn please create a Github issue including how we can replicate the issue and the expected vs actual behavior.


Contributing
============

Detailed instructions on developer setup and how to contribute are available `in the repo <https://github.com/bio-learn/biolearn/blob/master/DEVELOPMENT.md>`_
