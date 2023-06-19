# Bio-learn

## Table of Contents
- [Bio-learn](#bio-learn)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
  - [Pipenv shell](#pipenv-shell)
  - [Code Formatting with Black](#code-formatting-with-black)
  - [Testing with Pytest](#testing-with-pytest)
  - [Creating a Pull Request](#creating-a-pull-request)
  - [Working with Jupyter Notebook](#working-with-jupyter-notebook)

## Introduction

Bio-learn is a library for working with [biomarkers](https://en.wikipedia.org/wiki/Biomarker). Specifcally the library supports the ability to easily load data from disparate public sources into a common data structure, run predictive functions against that data and then plot the results for viewing.

## Getting Started

### Prerequisites

1. Python 3.8+
2. Pipenv
    ```
    pip install pipenv
    ```

### Installation

1. Clone the repo:

    ```
    git clone https://github.com/Biomarkers-of-Aging-Consortium/bio-learn
    cd bio-learn
    ```
1. Navigate to library directory
    ```
    cd boa-learn
    ```
2. Install project dependencies:

    ```
    pipenv install --dev
    ```

    This will create a virtual environment and install all the dependencies specified in the Pipfile.

## Pipenv shell

The pipenv shell allows you to run commands within the context of the dependencies described in the Pipfile. This helps ensure consistent behavior across environments. All of the below commands prefaced with `pipenv run` can be run without that preface while withing the shell. (ie `pytest` withing the shell vs  `pipenv run pytest` outside the shell)

```
pipenv shell
```

## Code Formatting with Black

We use Black for code formatting. Run it before committing to ensure your code follows our formatting guidelines:

```
pipenv run black .
```

## Testing with Pytest

To run tests, use Pytest:

```
pipenv run pytest
```

Make sure your code passes all the tests before pushing or creating a Pull Request.

## Creating a Pull Request

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

Please ensure your PR title is descriptive, and in the PR description, explain your changes, why you made them, and link any relevant issues.

## Working with Jupyter Notebook

This library is packaged with Jupyter notebooks that provide samples of how to use the library. Currently getting the notebooks to run correctly is slightly involved

1. Ensure you are in the library directory `boa-learn`
2. Start the shell

    ```
    pipenv shell
    ```
3. Once in the shell navigate down to the repo root
    ```
    cd ..
    ```
4. Launch jupyter
    ```
    jupyter notebook
    ```
5. You can now visit `localhost:8888` in your browser to start working on your notebooks. Notebooks are in the `notebooks` directory

