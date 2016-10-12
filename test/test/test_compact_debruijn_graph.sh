#!/bin/sh

DIR=test

rm -rf ${DIR}

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		mkdir -p ${DIR}/${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${DIR}/${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s} > ${DIR}/${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${DIR}/${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${DIR}/${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${DIR}/${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s} > ${DIR}/${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${DIR}/${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s} > ${DIR}/${s}${i}_t.log 2>&1
		
	done

done

echo "compare unitiqs to original"
diff ${DIR}/q1/u_chain.fasta ~/src/bliss/test/data/test.unitiqs.fasta


echo "COMPARE TO 1 proc, high mem, fastq"

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		echo "compare ${s}${i}"
		diff ${DIR}/q1 ${DIR}/${s}${i}
		
	done

done
