#!/bin/zsh
# mux_video.sh -- join the director's silent video and narration track into
# the upload file. make_conjugate_video.m does this itself when ffmpeg is
# installed; use this after installing ffmpeg (brew install ffmpeg), or to
# re-join after swapping audio without re-rendering.
#   ./mux_video.sh            # full quality
#   ./mux_video.sh _preview   # the preview render
set -e
cd "$(dirname "$0")/out"
TAG=${1:-}
ffmpeg -y -loglevel error -i "conjugate_points_silent$TAG.mp4" -i "narration$TAG.wav" \
       -c:v copy -c:a aac -b:a 192k -shortest "conjugate_points_video$TAG.mp4"
echo "wrote out/conjugate_points_video$TAG.mp4"
