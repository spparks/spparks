files="*.h *.cpp"
for i in $files
do
	grep -q "Copyright (2008) Sandia Corporation" ${i}
	if [[ $? -eq 1 ]]; then
		echo "${i} does not contain a copyright statement -- adding one in..."
		cat copyright ${i} > temp && mv temp ${i} #| sed -i '1d' ${i}
	else
		echo "Updating copyright statement on file ${i}"
                sed -i '6,9d; 5a\   Copyright(C) 1999-2025 National Technology & Engineering Solutions\
                of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with\
                NTESS, the U.S. Government retains certain rights in this software.' ${i}
	fi
done
