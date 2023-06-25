Biolearn
========

Biolearn enables easy and versatile analyses of biomarkers of aging data. It provides statistics and machine-learning tools, with instructive documentation and an open science approach.


.. warning::

    At the moment we do not have a stable release, only a **unstable developer version**.
    So you have to download the repo from github and install locally with **"pip install .[dev]"**


Important links
===============

- Official source code repo: https://github.com/bio-learn/biolearn/
- HTML documentation (**only development release for the moment**): https://bio-learn.github.io/

Install
=======

Latest release
--------------
**0. Install the right Python version**

Biolearn supports Python 3.10 and 3.11, please verify that you have the right Python 3 version before proceeding with the developer setup.

**1. Setup a virtual environment**

At the moment we do not have a stable release, only a unstable developer version.
This means that you should not install the package directly from pypi but always get the source code from github and do a local pip install.

In the following steps we show an installation path we have tested with venv and pip. Other package managers have not been tested and are not supported at the moment.

In linux as example create a virtual environment with venv:

.. code-block:: bash

    python3 -m venv bioenv
    source bioenv/bin/activate

Assuming you move to the local repository folder, you can install the package in the activated virtual environment:

.. code-block:: bash

    pip install -e .[dev]

Remember we do not support a pypi version at the moment, as in (pip install biolearn), the command will work but you will just get a dummy package.

Check installation
------------------

Try importing biolearn in a python / iPython session:

.. code-block:: python

    import biolearn

If no error is raised, you have installed biolearn correctly.

Discord server
==============

The biolearn team has a discord server to answer questions,
discuss feature requests, or have any biolearn related discussions.

Dependencies
============

The required dependencies to use the software are listed in the file `pyproject.toml <https://github.com/bio-learn/biolearn/blob/master/pyproject.toml>`_.

If you are using biolearn plotting functionalities or running the examples, matplotlib >= 3.3.0 is required.

Some plotting functions in biolearn support both matplotlib and plotly as plotting engines.

If you want to run the tests, you need pytest >= 6.0.0 and pytest-cov for coverage reporting.

Development
===========

Detailed instructions on developer setup and how to contribute are available `in the repo <https://github.com/bio-learn/biolearn/blob/master/DEVELOPMENT.md>`_