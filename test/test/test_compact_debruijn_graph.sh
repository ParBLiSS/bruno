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

		
		mkdir -p tmp/r${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > tmp/r${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/r${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s}        > tmp/r${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > tmp/r${s}${i}_t.log 2>&1
		
		
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
		
	mkdir -p tmp/tr1q${i}
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 1 -O tmp/tr1q${i}/s ~/src/bliss/test/data/test.star.unitiqs.fastq   > tmp/tr1q${i}_s.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 1 -O tmp/tr1q${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > tmp/tr1q${i}_s2.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 1 -O tmp/tr1q${i}/u ~/src/bliss/test/data/test.unitiqs.fastq        > tmp/tr1q${i}_u.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 1 -O tmp/tr1q${i}/t ~/src/bliss/test/data/test.short.unitiqs.fastq  > tmp/tr1q${i}_t.log 2>&1

	mkdir -p tmp/tr2q${i}
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 2 -O tmp/tr2q${i}/s ~/src/bliss/test/data/test.star.unitiqs.fastq   > tmp/tr2q${i}_s.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 2 -O tmp/tr2q${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > tmp/tr2q${i}_s2.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 2 -O tmp/tr2q${i}/u ~/src/bliss/test/data/test.unitiqs.fastq 		 > tmp/tr2q${i}_u.log 2>&1
	/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fastq -T -L 2 -O tmp/tr2q${i}/t ~/src/bliss/test/data/test.short.unitiqs.fastq  > tmp/tr2q${i}_t.log 2>&1

done

echo "compare unitiqs to original"
diff tmp/q1/u_chain.fasta ~/src/bliss/test/data/test.unitiqs.fasta


sed -n 0~2p tmp/q1/s2_chain.fasta | sort > s2_chains.sorted.txt
sed -n 0~2p tmp/q1/s_chain.fasta | sort > s_chains.sorted.txt
sed -n 0~2p tmp/q1/u_chain.fasta | sort > u_chains.sorted.txt
sed -n 0~2p tmp/q1/t_chain.fasta | sort > t_chains.sorted.txt

sed -n 0~2p tmp/tl2q1/s2_chain.fasta | sort > tl2_s2_chains.sorted.txt
sed -n 0~2p tmp/tl2q1/s_chain.fasta | sort  > tl2_s_chains.sorted.txt


echo "COMPARE TO 1 proc, high mem, fastq"
echo "COMPARE TO 1 proc, low mem, filtered, fastq"


for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		echo "compare ${s}${i}"
		diff tmp/q1 tmp/${s}${i}
		
		echo "compare l${s}${i}"
		diff -x "*.fasta" -x "*.debug" tmp/q1 tmp/l${s}${i}


		sed -n 0~2p tmp/l${s}${i}/s2_chain.fasta | sort > s2sorted.txt
		diff s2_chains.sorted.txt s2sorted.txt
		sort tmp/l${s}${i}/s2_compressed_chain.debug > s2sortedc.txt
		diff s2_chains.sorted.txt s2sortedc.txt

		sed -n 0~2p tmp/l${s}${i}/s_chain.fasta | sort > ssorted.txt
		diff s_chains.sorted.txt ssorted.txt
		sort tmp/l${s}${i}/s_compressed_chain.debug > ssortedc.txt
		diff s_chains.sorted.txt ssortedc.txt

		sed -n 0~2p tmp/l${s}${i}/u_chain.fasta | sort > usorted.txt
		diff u_chains.sorted.txt usorted.txt
		sort tmp/l${s}${i}/u_compressed_chain.debug > usortedc.txt
		diff u_chains.sorted.txt usortedc.txt

		sed -n 0~2p tmp/l${s}${i}/t_chain.fasta | sort > tsorted.txt
		diff t_chains.sorted.txt tsorted.txt
		sort tmp/l${s}${i}/t_compressed_chain.debug > tsortedc.txt
		diff t_chains.sorted.txt tsortedc.txt


		echo "compare r${s}${i}"
		diff -x "*.fasta" -x "*.debug" tmp/q1 tmp/r${s}${i}


		sed -n 0~2p tmp/r${s}${i}/s2_chain.fasta | sort > s2sorted.txt
		diff s2_chains.sorted.txt s2sorted.txt
		sort tmp/r${s}${i}/s2_compressed_chain.debug > s2sortedc.txt
		diff s2_chains.sorted.txt s2sortedc.txt

		sed -n 0~2p tmp/r${s}${i}/s_chain.fasta | sort > ssorted.txt
		diff s_chains.sorted.txt ssorted.txt
		sort tmp/r${s}${i}/s_compressed_chain.debug > ssortedc.txt
		diff s_chains.sorted.txt ssortedc.txt

		sed -n 0~2p tmp/r${s}${i}/u_chain.fasta | sort > usorted.txt
		diff u_chains.sorted.txt usorted.txt
		sort tmp/r${s}${i}/u_compressed_chain.debug > usortedc.txt
		diff u_chains.sorted.txt usortedc.txt

		sed -n 0~2p tmp/r${s}${i}/t_chain.fasta | sort > tsorted.txt
		diff t_chains.sorted.txt tsorted.txt
		sort tmp/r${s}${i}/t_compressed_chain.debug > tsortedc.txt
		diff t_chains.sorted.txt tsortedc.txt

		
	done

	echo "compare tl1q${i}"
	diff -x "*.valid" -x "*.fasta" -x "*.debug" tmp/q1 tmp/tl1q${i}


	sed -n 0~2p tmp/tl1q${i}/s2_chain.fasta | sort > s2sorted.txt
	diff s2_chains.sorted.txt s2sorted.txt
		sort tmp/tl1q${i}/s2_compressed_chain.debug > s2sortedc.txt
		diff s2_chains.sorted.txt s2sortedc.txt

	sed -n 0~2p tmp/tl1q${i}/s_chain.fasta | sort > ssorted.txt
	diff s_chains.sorted.txt ssorted.txt
		sort tmp/tl1q${i}/s_compressed_chain.debug > ssortedc.txt
		diff s_chains.sorted.txt ssortedc.txt

	sed -n 0~2p tmp/tl1q${i}/u_chain.fasta | sort > usorted.txt
	diff u_chains.sorted.txt usorted.txt
		sort tmp/tl1q${i}/u_compressed_chain.debug > usortedc.txt
		diff u_chains.sorted.txt usortedc.txt

	sed -n 0~2p tmp/tl1q${i}/t_chain.fasta | sort > tsorted.txt
	diff t_chains.sorted.txt tsorted.txt
		sort tmp/tl1q${i}/t_compressed_chain.debug > tsortedc.txt
		diff t_chains.sorted.txt tsortedc.txt


	echo "compare tl2q${i}"
	diff -x "*_chain.fasta" -x "*.debug" tmp/tl2q1 tmp/tl2q${i}

	sed -n 0~2p tmp/tl2q${i}/s2_chain.fasta | sort > s2sorted.txt
	diff tl2_s2_chains.sorted.txt s2sorted.txt
		sort tmp/tl2q${i}/s2_compressed_chain.debug > s2sortedc.txt
		diff tl2_s2_chains.sorted.txt s2sortedc.txt

	sed -n 0~2p tmp/tl2q${i}/s_chain.fasta | sort > ssorted.txt
	diff tl2_s_chains.sorted.txt ssorted.txt
		sort tmp/tl2q${i}/s_compressed_chain.debug > ssortedc.txt
		diff tl2_s_chains.sorted.txt ssortedc.txt



	echo "compare tr1q${i}"
	diff -x "*.valid" -x "*.fasta" -x "*.debug" tmp/q1 tmp/tr1q${i}


	sed -n 0~2p tmp/tr1q${i}/s2_chain.fasta | sort > s2sorted.txt
	diff s2_chains.sorted.txt s2sorted.txt
		sort tmp/tr1q${i}/s2_compressed_chain.debug > s2sortedc.txt
		diff s2_chains.sorted.txt s2sortedc.txt

	sed -n 0~2p tmp/tr1q${i}/s_chain.fasta | sort > ssorted.txt
	diff s_chains.sorted.txt ssorted.txt
		sort tmp/tr1q${i}/s_compressed_chain.debug > ssortedc.txt
		diff s_chains.sorted.txt ssortedc.txt

	sed -n 0~2p tmp/tr1q${i}/u_chain.fasta | sort > usorted.txt
	diff u_chains.sorted.txt usorted.txt
		sort tmp/tr1q${i}/u_compressed_chain.debug > usortedc.txt
		diff u_chains.sorted.txt usortedc.txt

	sed -n 0~2p tmp/tr1q${i}/t_chain.fasta | sort > tsorted.txt
	diff t_chains.sorted.txt tsorted.txt
		sort tmp/tr1q${i}/t_compressed_chain.debug > tsortedc.txt
		diff t_chains.sorted.txt tsortedc.txt


	echo "compare tr2q${i}"
	diff -x "*_chain.fasta" -x "*.debug" tmp/tl2q1 tmp/tr2q${i}

	sed -n 0~2p tmp/tr2q${i}/s2_chain.fasta | sort > s2sorted.txt
	diff tl2_s2_chains.sorted.txt s2sorted.txt
		sort tmp/tr2q${i}/s2_compressed_chain.debug > s2sortedc.txt
		diff tl2_s2_chains.sorted.txt s2sortedc.txt

	sed -n 0~2p tmp/tr2q${i}/s_chain.fasta | sort > ssorted.txt
	diff tl2_s_chains.sorted.txt ssorted.txt
		sort tmp/tr2q${i}/s_compressed_chain.debug > ssortedc.txt
		diff tl2_s_chains.sorted.txt ssortedc.txt




done
