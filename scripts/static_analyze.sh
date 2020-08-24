#!/bin/bash

# Assuming CodeChecker is installed and sourced
source ~/codechecker/venv/bin/activate
export PATH=$PATH:~/codechecker/build/CodeChecker/bin/

# fire up service
CodeChecker server&

# check the code
CodeChecker check \
    -b "./Allwclean; wmake libso src/rsr" \
    -i scripts/analyzer_skipfile \
    --checker-config \
    clang-tidy:cppcoreguidelines-special-member-functions.AllowMissingMoveFunctions=1\
    -o /tmp/static_results

# store check results
CodeChecker store /tmp/static_results \
    -n openrsr \
    --url localhost:8001/OpenRSR

# Generate html reports
CodeChecker parse -e html /tmp/results -o ./reports_html

rm -rf /tmp/static_results

firefox http://localhost:8001
