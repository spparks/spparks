#!/bin/bash

SPPARKS=$HOME/jaks.git/spparks.github/src/spk_spencer.gnu

# Array of layer intervals / thicknesses
declare -a Z=("0 25" "25 50" "50 75" "75 100")

# Array of "am build start" values
declare -a START=("25" "50" "75" "100")

# Array of Z values corresponding with bottom of build layer
declare -a LAYER_BOTTOM=("0" "25" "50" "75")

# Array of cartesian_layer commands
declare -a LAYERS=( \
  "1 start LL pass_id 1 thickness 25 offset -80.0 0.0" \
  "2 start UL pass_id 2 thickness 25 offset 0.0 80.0" \
  "3 start UR pass_id 1 thickness 25 offset 80.0 0.0" \
  "4 start LR pass_id 2 thickness 25 offset 0.0 -80.0" )

# Array of AM simulations windows
declare -a D=("0 25" "0 50" "0 75" "10 100")

NUM_LAYERS=${#Z[@]}

function initialize_stitchfile () {
   # This function runs a potts script on first layer 
   #   to create and initialize the stitch file.
   # Only run run this potts script on very first 
   #    layer to initialize stitch file
   L=0
   #  
   # Substitute region z layer information into 'in.init'
   layer=${Z[$L]}
   cat in.init | sed s"/Z/${layer}/" > in.potts_init
   
   #  
   # Run SPPARKS to initialize microstructure on layer
   SEED=$RANDOM
   $SPPARKS -var SEED $SEED < in.potts_init
}

# Initalize stitch file by running above function
initialize_stitchfile

# Loop of layers and simulate build microstructure
for (( L=0; L<${NUM_LAYERS}; L++ )); do 

  #  
  # Simulation am build on layer
  window=${D[$L]}
  layer=${LAYERS[$L]}
  layer_bottom=${LAYER_BOTTOM[$L]}
  cat in.am | sed s"/LAYER/${layer}/" \
	    | sed s"/WINDOW/${window}/" \
            | sed s"/BOTTOM/${layer_bottom}/" > in.am_layer

  #  
  # Run SPPARKS to simulate microstructure on each layer
  SEED=$RANDOM
  $SPPARKS -var SEED $SEED < in.am_layer

done
