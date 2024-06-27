SHELL := /bin/bash

.PHONY: clean test

flist = $(wildcard pf2barcode/figures/figure*.py)
allOutput = $(patsubst pf2barcode/figures/figure%.py, output/figure%.svg, $(flist))

all: $(allOutput)

output/figure%.svg: pf2barcode/figures/figure%.py
	@ mkdir -p ./output
	poetry run fbuild $*

test:
	poetry run pytest -s -x -v pf2barcode/tests/test_import.py

coverage.xml:
	poetry run pytest --cov=pf2barcode --cov-report=xml

clean:
	rm -rf output profile profile.svg
	rm -rf factor_cache

testprofile:
	poetry run python3 -m cProfile -o profile -m pytest -s -v -x
	gprof2dot -f pstats --node-thres=5.0 profile | dot -Tsvg -o profile.svg

mypy:
	poetry run mypy --install-types --non-interactive --ignore-missing-imports --check-untyped-defs pf2barcode
