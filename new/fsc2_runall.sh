#!/usr/bin/bash

for i in {1..100}
do 
	n=1	
	mkdir -p Model${i}_Combined
	while [ ${n} -le 100 ]
	do
		cp Model0_MSFS.obs Model${i}Rep${n}_MSFS.obs | cp /home/prm/Desktop/fsc_linux64/Models/Model${i}/Model${i}.est /home/prm/Desktop/fsc_linux64/Model${i}Rep${n}.est | cp /home/prm/Desktop/fsc_linux64/Models/Model${i}/Model${i}.tpl /home/prm/Desktop/fsc_linux64/Model${i}Rep${n}.tpl
		./fsc25221 -t Model${i}Rep${n}.tpl -n 100000 -N 100000 -m -u -e Model${i}Rep${n}.est -M 0.001 -l 5 -L 5 -c 0
		cat Model${i}Rep${n}/Model${i}Rep${n}.bestlhoods >> Model${i}_Combined/Model${i}_results.txt
		((n++))
	done
done
rm Model*Rep*.tpl | rm Model*Rep*.est | rm Model*Rep*.par | rm Model*Rep*_MSFS.obs
echo "Simulations completed"



