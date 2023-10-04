#!/bin/bash
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
for target in "bad_apple"; do
  if ! "$scriptDir/$target/build.sh" $@; then
    exit $?
  fi
done
