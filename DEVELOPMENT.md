# Biolearn Developer Guide

## Table of Contents
- [Biolearn Developer Guide](#biolearn-developer-guide)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Setup](#setup)
  - [Code Formatting with Black](#code-formatting-with-black)
  - [Running tests](#running-tests)
  - [Documentation](#documentation)
    - [Build and run locally](#build-and-run-locally)
    - [Updating deployed site](#updating-deployed-site)
    - [Adding examples](#adding-examples)
  - [Creating a Pull Request](#creating-a-pull-request)
  - [Working with Jupyter Notebooks](#working-with-jupyter-notebooks)

## Introduction

Biolearn enables easy and versatile analyses of biomarkers of aging data. It provides statistics and machine-learning tools, with instructive documentation and an open science approach.

## Getting Started

### Prerequisites

1. Python 3.10+

### Setup
Setup your vitual environment at the expected location. Ensure you are using a version of python that is 3.10 or higher as lower versions are not supported.

  ```
  python3 -m venv .venv
  ```

Then install all the dependencies using make
  ```
  make install
  ```

Dependencies are defined in the [pyproject.toml](/pyproject.toml) and new ones can be added. Ensure you run `make install` again after adding any new dependencies.


## Code Formatting with Black

We use Black for code formatting. Run it before committing to ensure your code follows our formatting guidelines:

  ```
  make format
  ```

## Running tests

  ```
  make test
  ```

Make sure your code passes all the tests before pushing or creating a Pull Request.

## Documentation

Biolearn uses web based documentation which is deployed [here](https://bio-learn.github.io/).

### Build and run locally

To build the doc site locally follow these steps from the repo root:

  ```
  source .venv/bin/activate
  cd doc
  make html
  ```
You can then navigate to the html pages locally starting at `doc/_build/html/index.html`

### Updating deployed site

To update the deployed docs copy the contents of `doc/_build/html` into the doc repo at https://github.com/bio-learn/bio-learn.github.io

### Adding examples

To provide example code for user to view and run add your scripts as a `.py` file in the `examples` folder. Subfolders of the `examples` folder are used to group examples into collections.

If you want the example to show on the homepage add it to `doc/index.rst` following the other examples in that file


## Creating a Pull Request

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

Please ensure your PR title is descriptive, and in the PR description, explain your changes, why you made them, and link any relevant issues.

## Working with Jupyter Notebooks

This library is packaged with examples that can be downloaded via the doc website as Jupyter notebooks. In order to run these notebooks with all the expected dependencies in the vitrual environment simply:

  ```
  make jupyter
  ```

