#!/usr/bin/env bash
# Monta o site em _site/ (ou na pasta dada como argumento).
# Requer pandoc. Uso:  bash site/build.sh [pasta_de_saida]
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-_site}"
rm -rf "$OUT"
mkdir -p "$OUT"

cp site/index.html site/style.css site/site.js site/favicon.svg "$OUT"/
cp -r docs/figs "$OUT"/figs
touch "$OUT"/.nojekyll

pagina() {   # pagina <fonte.md> <saida.html> <textos.yaml>
  pandoc "$1" -f gfm -t html5 --standalone \
    --template site/modelo-tutorial.html --lua-filter site/filtros.lua \
    --toc --toc-depth=3 --mathml \
    --metadata-file "$3" -o "$OUT/$2"
}

pagina docs/tutorial.md    tutorial.html    site/pt.yaml
pagina docs/tutorial_en.md tutorial-en.html site/en.yaml
echo "Site gerado em $OUT/"
