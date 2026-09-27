#!/bin/bash
# Build the linetools documentation locally.
#   Run from the docs/ directory.
set -e

rm -rf _build/html
rm -rf api

sphinx-build -W -b html . _build/html
echo "Docs written to docs/_build/html/index.html"
