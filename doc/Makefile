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
GH_PAGES_REPO = ../../bio-learn.github.io/

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

generate_dynamic_content:
	@echo "Generating dynamic content..."
	python generate.py

html:	sym_links generate_dynamic_content
	# These two lines make the build a bit more lengthy, and the
	# the embedding of images more robust
	rm -rf $(BUILDDIR)/html/_images
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	touch $(BUILDDIR)/html/.nojekyll
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

deploy: html
	# Assuming the docs repo is cloned in the same repo as biolearn is this will update the site
	rsync -av --delete --exclude .git/ _build/html/ $(GH_PAGES_REPO)
	cd $(GH_PAGES_REPO) && git add . && git commit -m "Automation Updated documentation" && git push