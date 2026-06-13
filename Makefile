.PHONY: all test figures values notebooks book clean

PY ?= python3

all: test values figures

test:
	$(PY) -m pytest -q

values:
	$(PY) scripts/make_values.py

figures:
	$(PY) scripts/make_figures.py

notebooks:
	jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb

book:
	@echo "LaTeX build is a WP5 deliverable; see textbook/ for the skeleton."
	cd textbook && (latexmk -lualatex main.tex || echo "(LaTeX toolchain not installed)")

clean:
	rm -rf figures/*.png figures/*.pdf .pytest_cache **/__pycache__
