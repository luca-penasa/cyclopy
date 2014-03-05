#!/bin/bash
cd doc
sphinx-apidoc -o apidoc ../cyclopy -s rst
make html
cd ..

URL="doc/_build/html/index.html"

# from stackoverflow q:3124556
[[ -x $BROWSER ]] && exec "$BROWSER" "$URL"
path=$(which xdg-open || which gnome-open) && exec "$path" "$URL"
echo "Can't find browser"


