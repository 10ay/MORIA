#!/bin/bash

for src in *.F; do
    [ -e "$src" ] || continue

    exe="${src%.F}.xgf"
    echo "Compiling $src -> $exe"
    gfortran-mp-14 -g "$src" -o "$exe"

    if [ $? -ne 0 ]; then
        echo "❌ Error compiling $src"
        exit 1
    fi
done

echo "Done compiling all .F files."
