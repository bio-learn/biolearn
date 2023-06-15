# simple makefile to simplify repetitive build env management tasks under posix
# CAREFUL: several commands are not compatible with Windows shell

PYTHON ?= python

all: clean test

clean-pyc:
	find . -name "*.pyc" | xargs rm -f
	find . -name "__pycache__" | xargs rm -rf

clean-build:
	rm -rf build

clean: clean-build clean-pyc

test-code:
	python -m pytest --pyargs biolearn --cov=biolearn

test-coverage:
	rm -rf coverage .coverage
	pytest --pyargs biolearn --showlocals --cov=biolearn --cov-report=html:coverage

test: test-code
