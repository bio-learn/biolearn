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

Biolearn supports Python 3.10, please verify that you have the right Python 3 version before proceeding with the developer setup.

**1. Setup a virtual environment**

At the moment we do not have a stable release, only a unstable developer version.
This means that you should not install the package directly from pypi but always get the source code from github and do a local pip install.

In the following steps we show an installation path we have tested with venv and pip. Other package managers have not been tested and are not supported at the moment.

In linux from inside the repo directory run the following to setup a virtual environment and install dependencies

.. code-block:: bash

    python3.10 -m venv .venv
    make install


Remember we do not support a pypi version at the moment, as in (pip install biolearn), the command will work but you will just get a dummy package.

Check installation
------------------

You can verify that your setup is working by running the tests

.. code-block:: bash

    make test

If no error is raised, you have installed biolearn correctly.

You can also run the following to start jupyter lab with the library available

.. code-block:: bash

    make jupyter

.. code-block:: python

    import biolearn

Discord server
==============

The biolearn team has a discord server to answer questions,
discuss feature requests, or have any biolearn related discussions.

Dependencies
============

The required dependencies to use the software are listed in the file `pyproject.toml <https://github.com/bio-learn/biolearn/blob/master/pyproject.toml>`_.


Development
===========

Detailed instructions on developer setup and how to contribute are available `in the repo <https://github.com/bio-learn/biolearn/blob/master/DEVELOPMENT.md>`_