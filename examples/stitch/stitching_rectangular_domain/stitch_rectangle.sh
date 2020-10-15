#!/bin/bash

# Required Input Parameters; 
# ONLY path to SPPARKS executable needs to be defined
SPPARKS=$HOME/jaks.git/spparks.cloned/src/spk_tutka.gnu
SPPARKS=$HOME/jaks.git/spparks.cloned/src/spk_spencer.gnu
# Edit this number as necessary; script will run in 10 to 30 
#  minutes depending upon how many procs are used; 
NUM_PROCS=16

function stitch_domain {

local F
local var
local BOTTOM Z1 THICKNESS
local Z0_HAZ DEPTH_HAZ MELT_DEPTH
local AM_PATH

while read LINE; do

  F=$(echo $LINE | cut -d ' ' -f1)
  #echo $F

  case "$F" in

    layer)
      var=$(echo $LINE | cut -d ' ' -f2)
      if [ "$var" == "window" ]; then
	BOTTOM=$(echo $LINE | cut -d ' ' -f3)
	Z1=$(echo $LINE | cut -d ' ' -f4)
	let THICKNESS=${Z1}-${BOTTOM}

      else
	echo "INPUT ERROR: expected 'layer window'"
        exit 1
      fi
      var=$(echo $LINE | cut -d ' ' -f5)
      if [ "$var" == "haz" ]; then
	Z0_HAZ=$(echo $LINE | cut -d ' ' -f6)
	let DEPTH_HAZ=${Z1}-${Z0_HAZ}
      else
	echo "INPUT ERROR: expected 'layer window'"
        exit 1
      fi
      var=$(echo $LINE | cut -d ' ' -f8)
      if [ "$var" == "melt_depth" ]; then
	MELT_DEPTH=$(echo $LINE | cut -d ' ' -f9)
      else
	echo "INPUT ERROR: expected 'melt_depth'"
        exit 1
      fi
      if [ ${BOTTOM} -lt ${Z0_HAZ} ]; then
	err="INPUT ERROR: elevation of layer bottom " 
	err+="cannot be less than depth_haz; check 'PathGen' data."
	echo $err
        exit 1
      fi
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
      AM_PATH="$(echo $LINE | cut -d ' ' -f 2,3,4,5,6,7)"

      var=$(echo $LINE | cut -d ' ' -f8)
      if [ "$var" != "block" ]; then
	echo "INPUT ERROR: expected 'block'"
        exit 1
      fi
      BLOCK="$(echo $LINE | cut -d ' ' -f 9,10,11,12)"
      LAYER_WINDOW="${BLOCK} ${BOTTOM} ${Z1}"
      WINDOW="${BLOCK} ${Z0_HAZ} ${Z1}"
      echo "LAYER_WINDOW=${LAYER_WINDOW}"
      echo "LAYER THICKNESS=${THICKNESS}"
      echo "WINDOW=${WINDOW}"
      echo "PATH=${AM_PATH}"
      echo "MELT_DEPTH=${MELT_DEPTH}"
      #####################################
      # Run potts model to initialize layer
      cat in.am | sed s"/MELT_DEPTH/${MELT_DEPTH}/" \
                | sed s"/DEPTH_HAZ/${DEPTH_HAZ}/" \
                | sed s"/PATH/${AM_PATH}/" \
                | sed s"/WINDOW/${WINDOW}/" \
                | sed s"/THICKNESS/${THICKNESS}/" > in.am_layer

      # Run SPPARKS to initialize microstructure on layer
      SEED=$RANDOM
      mpiexec -np $NUM_PROCS $SPPARKS -var SEED $SEED < in.am_layer
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

}

# Set pathgen data file
PATHDATA=stitch_rectangle.dat

# Run script
stitch_domain $PATHDATA

