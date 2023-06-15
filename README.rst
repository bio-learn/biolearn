Biolearn
========

Biolearn leverages the `scikit-learn <https://scikit-learn.org>`_ Python toolbox for biomarkers analysis

Important links
===============

- Official source code repo: https://github.com/bio-learn/biolearn/
- HTML documentation (stable release): https://bio-learn.github.io/

Install
=======

Latest release
--------------

**1. Setup a virtual environment**

We recommend that you install ``biolearn`` in a virtual Python environment.
In linux as example you can create one with venv:

.. code-block:: bash

    python3 -m venv bioenv
    source bioenv/bin/activate

After you can install the package in the activated virtual environment:

.. code-block:: bash

    pip install -U biolearn

Development version
-------------------

Please find all development setup instructions in the
`contribution guide <https://bio-learn.github.io/stable/development.html#setting-up-your-environment>`_.

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

For more information and ways to engage with the biolearn team see
`How to get help <https://bio-learn.github.io/stable/development.html#how-to-get-help>`_.

Dependencies
============

The required dependencies to use the software are listed in the file `pyproject.toml <https://github.com/bio-learn/biolearn/blob/main/pyproject.toml>`_.

If you are using biolearn plotting functionalities or running the examples, matplotlib >= 3.3.0 is required.

Some plotting functions in biolearn support both matplotlib and plotly as plotting engines.

If you want to run the tests, you need pytest >= 6.0.0 and pytest-cov for coverage reporting.

Development
===========

Detailed instructions on how to contribute are available at
http://bio-learn.github.io/stable/development.html
