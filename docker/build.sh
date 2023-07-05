#!/bin/bash

script_dir=$(readlink -f $0 | xargs dirname)
docker build -t detection-schema:latest -f "${script_dir}/Dockerfile" "${script_dir}/.."
