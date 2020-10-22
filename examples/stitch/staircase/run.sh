#!/bin/bash

SPPARKS=$HOME/jaks.git/spparks.cloned/src/spk_flamer.gnu
SPPARKS=$HOME/jaks.git/spparks.cloned/src/spk_spencer.gnu

# Layer thickness = 30 microns ~ 3 sites
let LAYER_THICKNESS=3
let STEP_HEIGHT=48
# 16=48/3
let RUNS_PER_STEP=16
# Height=480
# 10=480/48
let NUM_STEPS=10

let DEPTH_HAZ=5
let Z0=0
let TOP=480

# 
X0=0
X1=100
# With each step, domain in xy plane gets a bit smaller along y-axis
let DY=20
Y1=280

# Array of cartesian_layer commands
declare -a CARTESIAN_LAYERS=( \
  "1 start LL pass_id 1 thickness ${LAYER_THICKNESS} offset -80.0 0.0 serpentine 1" \
  "2 start UL pass_id 2 thickness ${LAYER_THICKNESS} offset 0.0 80.0 serpentine 1" \
  "3 start UR pass_id 1 thickness ${LAYER_THICKNESS} offset 80.0 0.0 serpentine 1" \
  "4 start LR pass_id 2 thickness ${LAYER_THICKNESS} offset 0.0 -80.0 serpentine 1" )

let START=1
let L=1

for (( s=1; s<=${NUM_STEPS}; s++ )); do

  let Y0=(${s}-1)*DY

  let Z1=(${s}-1)*STEP_HEIGHT+$LAYER_THICKNESS

  for (( r=1; r<=${RUNS_PER_STEP}; r++ )); do

    let LAYER_CMD=(${r})%4
    LAYER=${CARTESIAN_LAYERS[${LAYER_CMD}-1]}

    if [ "${START}" -eq 1 ]; then
      let Z0=0
      let START=0
    else
      let Z0=$Z1-$DEPTH_HAZ
    fi

    WINDOW="${X0} ${X1} ${Y0} ${Y1} ${Z0} ${Z1}"
    let BOTTOM=${Z1}-${LAYER_THICKNESS}
    #echo "r=${r}; LAYER=${LAYER}"
    #echo "WINDOW=${WINDOW}"

    # Run potts model to initialize first layer
    #if [ $L -eq 1 ]; then
    #   cat in.init | sed s"/WINDOW/${WINDOW}/" > in.potts_init
    #   
    #   #  
    #   # Run SPPARKS to initialize microstructure on layer
    #   SEED=$RANDOM
    #   $SPPARKS -var SEED $SEED < in.potts_init
    #fi

    # Create input file specialized for layer
    cat in.am | sed s"/LAYER/${LAYER}/" \
              | sed s"/DEPTH_HAZ/${DEPTH_HAZ}/" \
              | sed s"/WINDOW/${WINDOW}/" \
              | sed s"/BOTTOM/${BOTTOM}/" > in.am_layer

    # Run spparks
    SEED=$RANDOM
    mpiexec -np 16 $SPPARKS -var SEED $SEED < in.am_layer

    let Z1=$Z1+$LAYER_THICKNESS
    let L=$L+1

  done

done
