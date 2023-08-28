#!/bin/bash
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
for target in "cieutils" "linalg" "geo" "fem" "ciegl"; do
  "$scriptDir/$target/build.sh" $@
done

