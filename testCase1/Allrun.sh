#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

source $WM_PROJECT_DIR/bin/tools/RunFunctions

# create case
runApplication -s "fluid" blockMesh -region fluid
runApplication -s "fluid" snappyHexMesh -region fluid
runApplication -s "fluid" topoSet -region fluid
runApplication -s "fluid" extrudeToRegionMesh -region fluid
runApplication -s "film" topoSet -region film
runApplication decomposePar -allRegions
touch case.foam
paraFoam -touchAll

# run case
runParallel -a $(getApplication) || exit 2

