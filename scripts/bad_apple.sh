#!/bin/bash

set -e

args="$@"
shift $#

buildDirectory="$(dirname "$(dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )")")"
executable=$(find $buildDirectory -name bad_apple -type f | head -1)

trap "rm -rf bad_apple_img*.png" 0 2 3 15
ffmpeg -i $(yt-dlp -f 133 "https://www.youtube.com/watch?v=UkgK8eUdpAo" --get-url --rm-cache-dir) -vf format=gray "bad_apple_img%04d.png"

ls bad_apple_img*.png | $executable $args
