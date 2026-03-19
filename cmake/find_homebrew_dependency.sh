#!/usr/bin/env bash

packageVersions="$(brew search --formula "/$1(@[0-9]+)?/")"
for packageVersion in $(echo $packageVersions | tr ' ' '\n' | sort -r | tr '\n' ' '); do
    if brew list "$packageVersion" >/dev/null 2>&1; then
        foundPackage="$packageVersion"
        brew --prefix "$packageVersion" | tr -d '[:space:]'
        exit 0
    fi
done

echo "Error: no installed version of '$1' was found."
echo "Consider running 'brew install $1'."
exit 1
