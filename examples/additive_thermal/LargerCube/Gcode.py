# coding: utf-8
#This script converts slic3r generated g-code using the generic "RepRap" profile to a 3 column list of xyz points
#simulation heat source behavior will be determined by linearly interpolating between the specified points
#The simulation can then use this to determine what two points it is inbetween at any time

#The resulting path file is a 5-column file with lines of:
#X, Y, Z, Distance, PauseFlag.
#Pause flag values are: 0 (no pause), 1 (long, recoating pause) and 2 (beam movement pause)

import io
import numpy as np

fileName = "20mm_cube.gcode"
outputFile = "CubeEx.txt"
xMax = 300 #Use this value to scale the gcode simulation to the SPPARKS voxel domain.



#Actually reads in the file. Needs to disregard useless lines, strip off letter prefixes and keep track of Z at each line
def gCodeReader(fileName, mode = "r"):
	
	#I think we should actually use a regular python file reader, since all we want is every individual string from the file.
	with open(fileName, mode) as fIn:
		
		stringList = fIn.read().split()
				
	return stringList
		
		
#We have a rough setup for our input function, now we need to read in each string of the array, see what its first character is and then determine its meaning
def stringInterpreter(stringIn, stripChar):
	
	#we can call startswith to figure out the first character in the string. We can then test this against the specified stripChar (or x, y, z if we want to figure it out). Then well know whether to strip or not.
	if (stringIn.startswith(stripChar) == True):
		

		floatEquiv = float(stringIn.lstrip(stripChar))

		return floatEquiv
	#Our XY values should never be negative
	else:
		return -1

#This should loop through our input array, apply stringInterpreter to each value, figure out which are valid, and then put them in a new 5 column array. 
#The "F" values in the GCode file set the movement speed of the extruder. A value of F7800.000 appears to "max" out the speed. When this appears, we should
#turn on the pause flag. 
def arrayBuilder(inArray):

	zValue = 0
	zFlag = 0
	outArray = np.zeros((num_lines,5))
	k = 1 #Row counter for output array
	charTestLists = ['X','Y','Z','G','F']
	outArray[0,0:3] = 1 #Initialize first row with zeros
	outArray[0,3:4] = 0
	
	for i in range(0, len(inArray)):
		#print i
		#We need an inner loop to go through the possible values of stripChar.
		#Lets build an array of possible lead characters and then loop through the whole length.
		#We don't directly give the x,y,z assignment orders, but it should work out if things read in properly
		for j in charTestLists:
			
			if(len(inArray[i]) > 2):
				inValue = stringInterpreter(inArray[i],j)
			else:
				continue
				
			if( inValue == -1):
				
				continue
			else:
			
				if(j == 'X'):
					
					outArray[k,0] = inValue
					outArray[k,2] = zValue
					zFlag = 0
					
				elif(j == 'Y'):
				
					outArray[k,1] = inValue
					k = k + 1
					
				elif(j == 'Z'):
				
					zValue = np.floor(inValue)
					zFlag = 1
					
				#Let's not use G92 as our flag, as this isn't always called when we should have a pause
				#Instead, we'll look for F7800.000, which is when the nozzle is moving as fast as possible.
				#This will always be the last value on a line, and we should do the different types of pauses depending
				#on the values infront of it. If a F7800 line begins with Z, we should do a long pause and copy the X,Y values from the previous index
				#If it doesn't have a Z infront, we should assign the line a "short" pause.
				elif(j == 'F' and np.floor(inValue) == 7800):

					if(zFlag == 1):
						outArray[k,0] = outArray[k - 1, 0]
						outArray[k,1] = outArray[k - 1, 1]
						outArray[k,2] = zValue
						outArray[k - 1,4] = 1
						k = k + 1
					elif(outArray[k - 2,4] != 1):
						
						outArray[k - 2, 4] = 2
												
	#Now that we have all of our values, let's normalize things.
	#See where the array ends
	lastRow = 0
	for w in range(1,np.size(outArray[:,0])):
		if(outArray[w,0] == 0.0):
			lastRow = w
			break
	
	trimArray = outArray[2:lastRow,:]
	print trimArray
	#print lastRow
	#We need to use same normalization factor for all directions
	xDist = np.max(np.abs(trimArray[:,0])) - np.min(np.abs(trimArray[:,0]))
	print np.min(np.abs(trimArray[:,0])) 
# 	print np.min(trimArray[:,1])
	print xDist
	trimArray[:,0] = xMax/xDist * (np.abs(trimArray[:,0]) - np.min(np.abs(trimArray[:,0])))* (trimArray[:,0]/np.abs(trimArray[:,0]))
	trimArray[:,1] = xMax/xDist * (trimArray[:,1] - np.min(trimArray[:,1]))
	trimArray[:,2] = xMax/xDist * trimArray[:,2] + 5
	
	return trimArray
	
	
#We've filled in the x,y,z columns of our array, now lets calculate the cumulative distance traveled
#We do need to check if the machine ever stops and stops within a layer as this distance should not count!
def lengthCalculator(inArray):
	
	#Loop through the array and check if consecutive Z-values are the same. If they are, find the distance from X & Y.
	#If not, don't add any distance
	totalDistance = 0
	
	
	#Start at the second value, as the first will always be zero!
	for i in range(1,np.size(inArray,0)):
		
		if(inArray[i - 1,4] == 0):
			
			totalDistance = totalDistance + np.sqrt(pow(np.abs(inArray[i,0]) - np.abs(inArray[i-1,0]),2) + pow((inArray[i,1] - inArray[i-1,1]),2))
		
		inArray[i,3] = totalDistance
	
	return inArray
#Ok lets start testing the program

#Figure out how many lines are in the file
num_lines = sum(1 for line in open(fileName))

#Create a zero array of the correct size
finalArray = np.zeros((num_lines,4))

inputA = gCodeReader(fileName)
midArray = arrayBuilder(inputA)
#print(midArray[5000:5100,:])
finalArray = lengthCalculator(midArray)


#Output our final array to a new valeu
with open(outputFile, 'w') as f:
	print("XMax ", np.max(finalArray[:,0]))
	print('YMax', np.max(finalArray[:,1]))
	print('ZMax', np.max(finalArray[:,2]))
	np.savetxt(f, finalArray[:,:], fmt = '%10.3f', delimiter = ",")
	f.close()

