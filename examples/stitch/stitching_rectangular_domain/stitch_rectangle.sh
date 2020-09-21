#!/bin/bash

# Required Input Parameters
SPPARKS=$HOME/jaks.git/spparks.cloned/src/spk_spencer.gnu
SPPARKS=$HOME/jaks.git/spparks.cloned/src/spk_tutka.gnu
PATHDATA=stitch_rectangle.dat

declare -a THICKNESS
declare -a BOTTOM
declare -a Z1
declare -a Z0_HAZ
declare -a DEPTH_HAZ
declare -a MELT_DEPTH
declare -a LAYER_WINDOWS
declare -a LAYER_PTR

function get_layer_data {

let local L=0
let local PTR=0
local F
local var
local T

while read LINE; do

  F=$(echo $LINE | cut -d ' ' -f1)
  #echo $F

  case "$F" in

    layer)
      var=$(echo $LINE | cut -d ' ' -f2)
      if [ "$var" == "window" ]; then
	BOTTOM+=($(echo $LINE | cut -d ' ' -f3))
	Z1+=($(echo $LINE | cut -d ' ' -f4))
	let T=${Z1[$L]}-${BOTTOM[$L]}
	THICKNESS+=(${T})

      else
	echo "INPUT ERROR: expected 'layer window'"
        exit 1
      fi
      var=$(echo $LINE | cut -d ' ' -f5)
      if [ "$var" == "haz" ]; then
	Z0_HAZ+=($(echo $LINE | cut -d ' ' -f6))
	let T=${Z1[$L]}-${Z0_HAZ[$L]}
	DEPTH_HAZ+=(${T})
      else
	echo "INPUT ERROR: expected 'layer window'"
        exit 1
      fi
      var=$(echo $LINE | cut -d ' ' -f8)
      if [ "$var" == "melt_depth" ]; then
	MELT_DEPTH+=($(echo $LINE | cut -d ' ' -f9))
      else
	echo "INPUT ERROR: expected 'melt_depth'"
        exit 1
      fi
      if [ ${BOTTOM[$L]} -lt ${Z0_HAZ[$L]} ]; then
	err="INPUT ERROR: elevation of layer bottom " 
	err+="cannot be less than depth_haz; check 'PathGen' data."
	echo $err
        exit 1
      fi
      LAYER_PTR+=($PTR)
      let L+=1
      ;;
    path)
      # am path 1 start $X0 $Y0 end $X1 $Y1 speed 9
      # am path 1 PATH speed 9
      var=$(echo $LINE | cut -d ' ' -f2)
      if [ "$var" != "start" ]; then
	echo "INPUT ERROR: expected 'start'"
        exit 1
      fi
      var=$(echo $LINE | cut -d ' ' -f5)
      if [ "$var" != "end" ]; then
	echo "INPUT ERROR: expected 'end'"
        exit 1
      fi
      # Save this 'path' to array of PATHS
      PATHS+=("$(echo $LINE | cut -d ' ' -f 2,3,4,5,6,7)")

      var=$(echo $LINE | cut -d ' ' -f8)
      if [ "$var" != "block" ]; then
	echo "INPUT ERROR: expected 'block'"
        exit 1
      fi
      LAYER_WINDOWS+=("$(echo $LINE | cut -d ' ' -f 9,10,11,12)")
      # Now move PTR in preparation for next PATH
      let PTR+=1
      ;;
    \#)
         #echo COMMENT: ${LINE}
      ;;
    *)
      echo "Usage: {layer|domain}"
      exit 1
      ;;

  esac

done < $1

# Add one extra layer pointer to start of non existent next layer
LAYER_PTR+=($PTR)

}

function stitch {

let local NUM_LAYERS=${#THICKNESS[@]}

echo "NUM_LAYERS=${NUM_LAYERS}"
let local PTR=0
for (( L=0; L<$NUM_LAYERS; L++ )); do
  
  # Number of paths on layer
  let NUM_PATHS=${LAYER_PTR[$L+1]}-${LAYER_PTR[$L]}
  echo "Layer # ${LAYNUM}; Number of paths on layer = ${NUM_PATHS}"
  # Initialize layer
  for (( P=${LAYER_PTR[$L]}; P<${LAYER_PTR[$L+1]}; P++ )); do
    LAYER_WINDOW="${LAYER_WINDOWS[$P]} ${BOTTOM[$L]} ${Z1[$L]}"
    #echo "LAYER_WINDOW=${LAYER_WINDOW}"
    #####################################
    # Run potts model to initialize layer
    cat in.init | sed s"/WINDOW/${LAYER_WINDOW}/" > in.potts_init
    
    #  
    # Run SPPARKS to initialize microstructure on layer
    SEED=$RANDOM
    $SPPARKS -var SEED $SEED < in.potts_init
    #####################################
  done

  # Run AM model on layer
  for (( P=${LAYER_PTR[$L]}; P<${LAYER_PTR[$L+1]}; P++ )); do
    WINDOW="${LAYER_WINDOWS[$P]} ${Z0_HAZ[$L]} ${Z1[$L]}"
    AM_PATH=${PATHS[$P]}
    echo "WINDOW=${WINDOW}"
    echo "PATH=${AM_PATH}"
    echo "MELT_DEPTH=${MELT_DEPTH[$L]}"
    #####################################
    # Run potts model to initialize layer
    cat in.am | sed s"/MELT_DEPTH/${MELT_DEPTH[$L]}/" \
              | sed s"/DEPTH_HAZ/${DEPTH_HAZ[$L]}/" \
              | sed s"/PATH/${AM_PATH}/" \
              | sed s"/WINDOW/${WINDOW}/" \
              | sed s"/THICKNESS/${THICKNESS[$L]}/" > in.am_layer

    # Run SPPARKS to initialize microstructure on layer
    SEED=$RANDOM
    $SPPARKS -var SEED $SEED < in.am_layer
    #####################################
  done

done
}

get_layer_data $PATHDATA
stitch

