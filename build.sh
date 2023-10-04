#!/bin/bash
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
for directory in "libraries" "executables"; do
  if ! "$scriptDir/$directory/build.sh" $@; then
    exit 1
  fi
done

