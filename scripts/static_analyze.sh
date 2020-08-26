#!/bin/bash

# Assuming CodeChecker is installed and sourced

branch="$1"
commit="$2"
if [ -z "$1" ]; then
    branch="develop"
    commit="0000000"
fi

source /opt/openfoam7/etc/bashrc

set -x 

# check the code
CodeChecker check \
    -b "./Allwclean; ./Allwmake" \
    -i scripts/analyzer_skipfile \
    --checker-config \
    clang-tidy:cppcoreguidelines-special-member-functions.AllowMissingMoveFunctions=1\
    -o /tmp/static_results

# store check results
CodeChecker store /tmp/static_results \
    -n $commit --tag $branch \
    --url http://openrsr-code-check.herokuapp.com:80/OpenRSR
