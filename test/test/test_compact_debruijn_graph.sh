#!/bin/sh

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		mkdir -p tmp/${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s} > tmp/${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s} > tmp/${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s} > tmp/${s}${i}_t.log 2>&1
		
		mkdir -p tmp/l${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s} > tmp/l${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/l${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s} > tmp/l${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s} > tmp/l${s}${i}_t.log 2>&1

		
		
		
	done

	mkdir -p tmp/tl1q${i}
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1q${i}/s ~/src/bliss/test/data/test.star.unitiqs.fastq > tmp/tl1q${i}_s.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1q${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > tmp/tl1q${i}_s2.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1q${i}/u ~/src/bliss/test/data/test.unitiqs.fastq > tmp/tl1q${i}_u.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1q${i}/t ~/src/bliss/test/data/test.short.unitiqs.fastq > tmp/tl1q${i}_t.log 2>&1

	mkdir -p tmp/tl2q${i}
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q${i}/s ~/src/bliss/test/data/test.star.unitiqs.fastq > tmp/tl2q${i}_s.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > tmp/tl2q${i}_s2.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q${i}/u ~/src/bliss/test/data/test.unitiqs.fastq > tmp/tl2q${i}_u.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q${i}/t ~/src/bliss/test/data/test.short.unitiqs.fastq > tmp/tl2q${i}_t.log 2>&1
		

done

echo "compare unitiqs to original"
diff tmp/q1/u_chain.fasta ~/src/bliss/test/data/test.unitiqs.fasta


echo "COMPARE TO 1 proc, high mem, fastq"
echo "COMPARE TO 1 proc, low mem, filtered, fastq"

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		echo "compare ${s}${i}"
		diff tmp/q1 tmp/${s}${i}
		
		echo "compare l${s}${i}"
		diff -x "*_branch.fasta" tmp/q1 tmp/l${s}${i}
		
	done

	echo "compare tl1q${i}"
	diff -x "*.valid" -x "*_branch.fasta" tmp/q1 tmp/tl1q${i}

	echo "compare tl2q${i}"
	diff tmp/tl2q1 tmp/tl2q${i}


done
