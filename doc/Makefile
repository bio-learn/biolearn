# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = -v
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"

clean:
	-rm -rf $(BUILDDIR)/*
	-rm -rf auto_examples/
	-rm -rf generated/*
	-rm -rf modules/generated/*
	# -rm -rf ../examples/**/biolearn_cache
	# -rm -rf ../**/biolearn_cache


sym_links:
	# Make sym-links to share the cache across various example
	# directories
	# The following find command selects all the depth-one
	# subdirectories of ../examples/ and creates a sym-link to
	# ../examples/biolearn_cache/
	# -find ../examples/ -type d -mindepth 1 -maxdepth 1 -not -path "../examples/biolearn_cache" -exec ln -sf ../biolearn_cache "{}/biolearn_cache" \;
	# -mkdir -p ../examples/biolearn_cache/joblib

force_html: force html

force:
	find . -name \*.rst -exec touch {} \;

html:	sym_links
	# These two lines make the build a bit more lengthy, and the
	# the embedding of images more robust
	rm -rf $(BUILDDIR)/html/_images
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	touch $(BUILDDIR)/html/.nojekyll
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
