#!/bin/bash

PATHDATA=stitching.dat
declare -a PATHS
declare -a WINDOWS
declare -a LAYER_WINDOWS
let ZSTART=0
let DEPTH_HAZ=5


function get_paths_and_windows {

local F
local var
local Z0
let Z1=$ZSTART
local LAYER_BOTTOM

let local L=1
while read LINE; do

  F=$(echo $LINE | cut -d ' ' -f1)
  echo $F

  case "$F" in

    layer)
      var=$(echo $LINE | cut -d ' ' -f2)
      if [ "$var" == "thickness" ]; then
        let THICKNESS=$(echo $LINE | cut -d ' ' -f3)
      else
	echo "INPUT ERROR: expected 'layer thickness'"
        exit 1
      fi
      let Z1=$Z1+$THICKNESS
      if [ $L -ge 2 ]; then
	let Z0=$Z1-$DEPTH_HAZ
      else
	let Z0=$ZSTART
      fi
      let L=$L+1
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
      let LAYER_BOTTOM=${Z1}-${THICKNESS}
      
      LAYER_WINDOWS+=("$(echo $LINE | cut -d ' ' -f 9,10,11,12) ${LAYER_BOTTOM} ${Z1}")
      WINDOWS+=("$(echo $LINE | cut -d ' ' -f 9,10,11,12) ${Z0} ${Z1}")
      ;;
      \#)
         echo COMMENT: ${LINE}
      ;;
    *)
      echo "Usage: {layer|domain}"
      exit 1
      ;;

  esac

done < $1

}


get_paths_and_windows $PATHDATA

NUM_RUNS=${#PATHS[@]}

for (( r=0; r<$NUM_RUNS; r++ )); do
  echo LAYER_WINDOW=${LAYER_WINDOWS[$r]}
  echo WINDOW=${WINDOWS[$r]}
  echo PATH=${PATHS[$r]}

    # Run potts model to initialize first layer
    if [ $L -eq 1 ]; then
       cat in.init | sed s"/WINDOW/${WINDOW}/" > in.potts_init
       
       #  
       # Run SPPARKS to initialize microstructure on layer
       SEED=$RANDOM
       $SPPARKS -var SEED $SEED < in.potts_init
    fi
done
