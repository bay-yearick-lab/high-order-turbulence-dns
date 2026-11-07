# Python reference implementation: verification, tests, benchmarks, figures.
# The Fortran channel solver has its own makefile in fortran/.
PY ?= python

.PHONY: verify test bench assets

verify:          ## headline numbers in about 30 s on a laptop
	PYTHONPATH=python $(PY) -m bspline_filter.verify

test:            ## numerical validation suite
	$(PY) -m pytest -q tests

bench:           ## timing and iteration benchmarks
	$(PY) benchmarks/run_benchmarks.py

assets:          ## regenerate README figures in assets/
	$(PY) tools/make_readme_assets.py
