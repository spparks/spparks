test:
	rm -f *.st *.st-journal
	PYTHONPATH=..:../../../local/lib:${PYTHONPATH} python3 -X faulthandler unit_cv_readwrite.py
#	PYTHONPATH=..:../../../local/lib:${PYTHONPATH} mpirun -np 4 python3 -X faulthandler unit_cv_mp_readwrite.py 2 2 1
#	PYTHONPATH=..:../../../local/lib:${PYTHONPATH} mpirun -np 4 python3 -X faulthandler unit_potts_weld.py

clean:
	rm -f *.st *.st-journal
