#!/bin/bash
# RUN THIS SCRIPT FROM CASE DIR
# $1 Is supposed to be solver name

ext="png"


pwd
for ca in {1..1}; do
    sol=\'"$1"\'
    gnuplot -e "solver=$sol" -e "case=$ca" -e "start=0.0" plots/plotSaturationProfile.gp
    gnuplot -e "solver=$sol" -e "case=$ca" -e "start=0.5" plots/plotSaturationProfile.gp
    gnuplot -e "solver=$sol" -e "case=$ca" -e "start=0.0" plots/plotPressureProfile.gp
    gnuplot -e "solver=$sol" -e "case=$ca" -e "start=0.5" plots/plotPressureProfile.gp
    gnuplot -e "solver=$sol" -e "case=$ca" plots/plotSaturationAtInterface.gp
done
