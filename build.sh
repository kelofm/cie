#!/bin/bash
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
for directory in "libraries" "executables"; do
  "$scriptDir/$directory/build.sh" $@
done

