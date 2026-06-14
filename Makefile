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
	$(PY) scripts/make_values_tex.py

params:
	$(PY) scripts/make_param_study.py

figures:
	$(PY) scripts/make_figures.py

nbgen:
	$(PY) scripts/make_notebooks.py
	$(PY) scripts/make_textbook_skeleton.py

notebooks:
	jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb

book: figures
	$(PY) scripts/make_values_tex.py
	cd textbook && lualatex -interaction=nonstopmode main.tex >/dev/null \
	  && lualatex -interaction=nonstopmode main.tex >/dev/null \
	  && echo "built textbook/main.pdf"

clean:
	rm -rf figures/*.png figures/*.pdf .pytest_cache **/__pycache__
