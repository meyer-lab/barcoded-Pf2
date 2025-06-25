SHELL := /bin/bash

.PHONY: clean

flist = $(wildcard pf2barcode/figures/figure*.py)
allOutput = $(patsubst pf2barcode/figures/figure%.py, output/figure%.svg, $(flist))

all: $(allOutput)

.venv: pyproject.toml
	rye sync

clean:
	rm -r quarto/_output

pyright: .venv
	rye run pyright pf2barcode
