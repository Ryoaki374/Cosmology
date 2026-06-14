.PHONY: all test figures values params polarization class nbgen notebooks book clean

PY ?= python3

all: test values params figures nbgen

polarization:
	$(PY) scripts/make_polarization.py

class:
	$(PY) scripts/make_class_comparison.py

test:
	$(PY) -m pytest -q

values:
	$(PY) scripts/make_values.py

params:
	$(PY) scripts/make_param_study.py

figures:
	$(PY) scripts/make_figures.py

nbgen:
	$(PY) scripts/make_notebooks.py
	$(PY) scripts/make_textbook_skeleton.py

notebooks:
	jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb

book:
	@echo "LaTeX build is a WP5 deliverable; see textbook/ for the skeleton."
	cd textbook && (latexmk -lualatex main.tex || echo "(LaTeX toolchain not installed)")

clean:
	rm -rf figures/*.png figures/*.pdf .pytest_cache **/__pycache__
