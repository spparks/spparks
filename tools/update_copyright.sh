files="*.h *.cpp"
for i in $files
do
	line1=`grep -n "\/\*" ${i} | awk '{print($1)}' | sed 's/[^0-9]//g' | head -n1`
	line2=`grep -n "\*\/" ${i} | awk '{print($1)}' | sed 's/[^0-9]//g' | head -n1`
	ackLines=`grep -n Sandia ${i} | awk '{print($1)}' | sed 's/[^0-9]//g'`

	# Initialize a flag to true
	all_between=true

	# Loop through each value in lines
	for line in $ackLines; do
	    # Check if the acknowledgement statement is within expected bounds
	    if (( line < line1 || line > line2 )); then
		all_between=false
		break  # Exit the loop if any value is out of range
	    fi
	done

#	Uncomment to check for non-conforming files
#	if $all_between; then
#	    echo "All values in lines are between $line1 and $line2."
#	else
#	    echo "Not all values in lines are between $line1 and $line2 for file ${i}."
#	fi

	echo "Updating copyright statement on file ${i}"
	sed -i "${line1},${line2}d" ${i}
	cat copyright ${i} > temp && mv temp ${i}
done
