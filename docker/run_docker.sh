#!/usr/bin/env bash

indir=/path/to/input/dir
outdir=/path/to/output/dir

docker run --rm -i \
    -v "${indir}":/input \
    -v "${outdir}":/output \
    detection-schema:latest \
    -d /input \
    -o /output \
    -i
