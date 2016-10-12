#!/bin/sh

rm -rf tmp

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		mkdir -p tmp/${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s} > tmp/${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s} > tmp/${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s} > tmp/${s}${i}_t.log 2>&1
		
	done

done

echo "compare unitiqs to original"
diff tmp/q1/u_chain.fasta ~/src/bliss/test/data/test.unitiqs.fasta


echo "COMPARE TO 1 proc, high mem, fastq"

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		echo "compare ${s}${i}"
		diff tmp/q1 tmp/${s}${i}
		
	done

done
