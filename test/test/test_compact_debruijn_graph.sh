#!/bin/sh

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		mkdir ${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${s}${i}/s /home/tpan/src/bliss/test/data/test.star.unitiqs.fast${s} > ${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${s}${i}/s2 /home/tpan/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${s}${i}/u /home/tpan/src/bliss/test/data/test.unitiqs.fast${s} > ${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O ${s}${i}/t /home/tpan/src/bliss/test/data/test.short.unitiqs.fast${s} > ${s}${i}_t.log 2>&1
		
		mkdir l${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O l${s}${i}/s /home/tpan/src/bliss/test/data/test.star.unitiqs.fast${s} > l${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O l${s}${i}/s2 /home/tpan/src/bliss/test/data/test.star2.unitiqs.fast${s} > l${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O l${s}${i}/u /home/tpan/src/bliss/test/data/test.unitiqs.fast${s} > l${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O l${s}${i}/t /home/tpan/src/bliss/test/data/test.short.unitiqs.fast${s} > l${s}${i}_t.log 2>&1
		
	done
done

echo "COMPARE TO 1 proc, high mem, fastq"

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
	
		diff q1 ${s}${i}
		
		diff q1 l${s}${i}
		
	done
done
