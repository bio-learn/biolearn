# simple makefile to simplify repetitive build env management tasks under posix
# CAREFUL: several commands are not compatible with Windows shell
VENV_BIN     = .venv/bin
PYTHON   = $(VENV_BIN)/python
PIP      = $(VENV_BIN)/pip
PYTEST   = $(VENV_BIN)/pytest
BLACK    = $(VENV_BIN)/black
JUPYTER  = $(VENV_BIN)/jupyter
HATCH  = $(VENV_BIN)/hatch
FILENAME ?=

all: clean test format

clean-pyc:
	find . -name "*.pyc" | xargs rm -f
	find . -name "__pycache__" | xargs rm -rf

clean-build:
	rm -rf build

clean: clean-build clean-pyc

generate-test-data:
	$(PYTHON) ./biolearn/test/generate.py

jupyter: 
	mkdir -p notebooks
	$(JUPYTER) lab

install: 
	$(PIP) install -e .[dev]

format:
	$(BLACK) ./biolearn

check-format:
	$(BLACK) --check ./biolearn

publish: 
	$(HATCH) build -c
	$(HATCH) publish

test-code:
	$(PYTEST) --pyargs biolearn --cov=biolearn $(if $(strip $(FILENAME)), -k $(FILENAME))

test-coverage:
	rm -rf coverage .coverage
	$(PYTEST) --pyargs biolearn --showlocals --cov=biolearn --cov-report=html:coverage

test: generate-test-data test-code
