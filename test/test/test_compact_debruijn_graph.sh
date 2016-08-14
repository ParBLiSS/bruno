#!/bin/sh

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		mkdir -p tmp/${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/s /home/tpan/src/bliss/test/data/test.star.unitiqs.fast${s} > tmp/${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/s2 /home/tpan/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/u /home/tpan/src/bliss/test/data/test.unitiqs.fast${s} > tmp/${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s} -O tmp/${s}${i}/t /home/tpan/src/bliss/test/data/test.short.unitiqs.fast${s} > tmp/${s}${i}_t.log 2>&1
		
		mkdir -p tmp/l${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/s /home/tpan/src/bliss/test/data/test.star.unitiqs.fast${s} > tmp/l${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/s2 /home/tpan/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/l${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/u /home/tpan/src/bliss/test/data/test.unitiqs.fast${s} > tmp/l${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fast${s} -O tmp/l${s}${i}/t /home/tpan/src/bliss/test/data/test.short.unitiqs.fast${s} > tmp/l${s}${i}_t.log 2>&1

		
		mkdir -p tmp/tl1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1${s}${i}/s /home/tpan/src/bliss/test/data/test.star.unitiqs.fastq > tmp/tl1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1${s}${i}/s2 /home/tpan/src/bliss/test/data/test.star2.unitiqs.fastq > tmp/tl1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1${s}${i}/u /home/tpan/src/bliss/test/data/test.unitiqs.fastq > tmp/tl1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 1 -O tmp/tl1${s}${i}/t /home/tpan/src/bliss/test/data/test.short.unitiqs.fastq > tmp/tl1${s}${i}_t.log 2>&1
	
		mkdir -p tmp/tl2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2${s}${i}/s /home/tpan/src/bliss/test/data/test.star.unitiqs.fastq > tmp/tl2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2${s}${i}/s2 /home/tpan/src/bliss/test/data/test.star2.unitiqs.fastq > tmp/tl2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2${s}${i}/u /home/tpan/src/bliss/test/data/test.unitiqs.fastq > tmp/tl2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2${s}${i}/t /home/tpan/src/bliss/test/data/test.short.unitiqs.fastq > tmp/tl2${s}${i}_t.log 2>&1
		
		
		
	done


done

echo "compare unitiqs to original"
diff tmp/q1/u_chain.fasta /home/tpan/src/bliss/test/data/test.unitiqs.fasta


echo "COMPARE TO 1 proc, high mem, fastq"
echo "COMPARE TO 1 proc, low mem, filtered, fastq"

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		echo "compare ${s}${i}"
		diff tmp/q1 tmp/${s}${i}
		
		echo "compare l${s}${i}"
		diff tmp/q1 tmp/l${s}${i}

		echo "compare tl1${s}${i}"
		diff tmp/q1 tmp/tl1${s}${i}
		
		echo "compare tl2${s}${i}"
		diff tmp/tl2q1 tmp/tl2${s}${i}
	done
	

done