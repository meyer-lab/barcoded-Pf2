SHELL := /bin/bash

.PHONY: clean test

flist = $(wildcard pf2barcode/figures/figure*.py)
allOutput = $(patsubst pf2barcode/figures/figure%.py, output/figure%.svg, $(flist))

all: $(allOutput)

output/figure%.svg: pf2barcode/figures/figure%.py
	@ mkdir -p ./output
	rye run fbuild $*

test:
	rye run pytest -s -x -v pf2barcode/tests/test_import.py

clean:
	rm -rf output profile profile.svg
	rm -rf factor_cache

pyright: .venv
	rye run pyright pf2

coverage.xml: .venv
	rye run pytest --cov=pf2barcode --cov-report xml:coverage.xml

testprofile:
	rye run python3 -m cProfile -o profile -m pytest -s -v -x
	gprof2dot -f pstats --node-thres=5.0 profile | dot -Tsvg -o profile.svg
