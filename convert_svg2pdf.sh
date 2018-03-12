#!/bin/bash
for f in *.svg; do inkscape $f --export-pdf=${f/%.svg/.pdf}; done

#pdftk *.pdf cat output RAG_TLX3_agg.pdf

