#!/bin/bash
# Build and Test the toolkit in a container based on OF-enabled docker image
# Assumes OpenRSR git repo is at $USER/OpenRSR-OF7

# First parameter to this script should point to the base docker image
base_image=$1
if [ -z $base_image ]; then
    base_image="foamscience/bionic-openfoam7"
fi

# Second parameter should point to the path of openFOAM's bashrc inside
# the image
bashrc_path=$2
if [ -z $bashrc_path ]; then
    bashrc_path="/opt/openfoam7/etc/bashrc"
fi


check_errs()
{
  # Parameter 1 is the return code
  # Parameter 2 is text to display on failure.
  if [ "${1}" -ne "0" ]; then
    echo "ERROR # ${1} : ${2}"
    exit "${1}"
  fi
}

# Source OpenFOAM
source $bashrc_path
set -ev

# Get to the directory, compile libraries (Opt mode)
cd $USER/OpenRSR-OF7; ./scripts/installPetsc.sh &&  ./Allwclean; ./Allwmake

# Run tests
echo "Testing Libraries:"
echo "------------------------------------------\n"
./Alltest
