#!/bin/bash
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
for target in "bad_apple"; do
  "$scriptDir/$target/build.sh" $@
done

