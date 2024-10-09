#!/bin/bash

scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$scriptDir"

TAGFILES=

for lib in cieutils linalg geo fem ciegl; do
    cd "libraries/$lib/docs"
    if ! ( cat doxyfile ; echo "TAGFILES=$TAGFILES" ) | doxygen -; then
        exit 1
    fi
    TAGFILES="$TAGFILES ../../$lib/docs/${lib}_doxygen_tagfile=../../../$lib/docs/html"
    cd "$scriptDir"
done
