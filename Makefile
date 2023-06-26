# simple makefile to simplify repetitive build env management tasks under posix
# CAREFUL: several commands are not compatible with Windows shell
VENV_BIN     = .venv/bin
PYTHON   = $(VENV_BIN)/python
PIP      = $(VENV_BIN)/pip
PYTEST   = $(VENV_BIN)/pytest
BLACK    = $(VENV_BIN)/black
JUPYTER  = $(VENV_BIN)/jupyter

all: clean test

clean-pyc:
	find . -name "*.pyc" | xargs rm -f
	find . -name "__pycache__" | xargs rm -rf

clean-build:
	rm -rf build

clean: clean-build clean-pyc

jupyter: 
	$(JUPYTER) notebook

setup: 
	python3 -m venv .venv
	$(PIP) install -e .[dev]

format:
	$(BLACK) ./biolearn

test-code:
	$(PYTEST) --pyargs biolearn --cov=biolearn

test-coverage:
	rm -rf coverage .coverage
	$(PYTEST) --pyargs biolearn --showlocals --cov=biolearn --cov-report=html:coverage

test: test-code
