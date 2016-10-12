#!/bin/sh

rm -rf tmp


mkdir -p tmp/q1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O tmp/q1/s ~/src/bliss/test/data/test.star.unitiqs.fastq >   tmp/q1_s.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O tmp/q1/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > tmp/q1_s2.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O tmp/q1/u ~/src/bliss/test/data/test.unitiqs.fastq > 		tmp/q1_u.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O tmp/q1/t ~/src/bliss/test/data/test.short.unitiqs.fastq >  tmp/q1_t.log 2>&1

mkdir -p tmp/tl2q1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q1/s ~/src/bliss/test/data/test.star.unitiqs.fastq > 	tmp/tl2q1_s.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q1/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > 	tmp/tl2q1_s2.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q1/u ~/src/bliss/test/data/test.unitiqs.fastq > 		tmp/tl2q1_u.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O tmp/tl2q1/t ~/src/bliss/test/data/test.short.unitiqs.fastq > 	tmp/tl2q1_t.log 2>&1


for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		
		mkdir -p tmp/r${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > tmp/r${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/r${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s}        > tmp/r${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -O tmp/r${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > tmp/r${s}${i}_t.log 2>&1
		

		mkdir -p tmp/tr1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 1 -O tmp/tr1${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > tmp/tr1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 1 -O tmp/tr1${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/tr1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 1 -O tmp/tr1${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s}        > tmp/tr1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 1 -O tmp/tr1${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > tmp/tr1${s}${i}_t.log 2>&1
	
		mkdir -p tmp/tr2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 2 -O tmp/tr2${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > tmp/tr2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 2 -O tmp/tr2${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > tmp/tr2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 2 -O tmp/tr2${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s} 		  > tmp/tr2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -T -L 2 -O tmp/tr2${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > tmp/tr2${s}${i}_t.log 2>&1

		
	done

done

echo "compare unitiqs to original"
diff tmp/q1/u_chain.fasta ~/src/bliss/test/data/test.unitiqs.fasta


sed -n 0~2p tmp/q1/s2_chain.fasta | sort > s2_chains.sorted.txt
sed -n 0~2p tmp/q1/s_chain.fasta | sort > s_chains.sorted.txt
sed -n 0~2p tmp/q1/u_chain.fasta | sort > u_chains.sorted.txt
sed -n 0~2p tmp/q1/t_chain.fasta | sort > t_chains.sorted.txt

sed -n 0~2p tmp/tl2q1/s2_chain.fasta | sort > tl2_s2_chains.sorted.txt
sed -n 0~2p tmp/tl2q1/s_chain.fasta | sort  > tl2_s_chains.sorted.txt

echo "COMPARE TO 1 proc, refactored, filtered, fastq"

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 

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

		echo "compare tr1${s}${i}"
		diff -x "*.valid" -x "*.fasta" -x "*.debug" tmp/q1 tmp/tr1${s}${i}
	
		sed -n 0~2p tmp/tr1${s}${i}/s2_chain.fasta | sort > s2sorted.txt
		diff s2_chains.sorted.txt s2sorted.txt
			sort tmp/tr1${s}${i}/s2_compressed_chain.debug > s2sortedc.txt
			diff s2_chains.sorted.txt s2sortedc.txt
	
		sed -n 0~2p tmp/tr1${s}${i}/s_chain.fasta | sort > ssorted.txt
		diff s_chains.sorted.txt ssorted.txt
			sort tmp/tr1${s}${i}/s_compressed_chain.debug > ssortedc.txt
			diff s_chains.sorted.txt ssortedc.txt
	
		sed -n 0~2p tmp/tr1${s}${i}/u_chain.fasta | sort > usorted.txt
		diff u_chains.sorted.txt usorted.txt
			sort tmp/tr1${s}${i}/u_compressed_chain.debug > usortedc.txt
			diff u_chains.sorted.txt usortedc.txt
	
		sed -n 0~2p tmp/tr1${s}${i}/t_chain.fasta | sort > tsorted.txt
		diff t_chains.sorted.txt tsorted.txt
			sort tmp/tr1${s}${i}/t_compressed_chain.debug > tsortedc.txt
			diff t_chains.sorted.txt tsortedc.txt

		echo "compare tr2${s}${i}"
		diff -x "*_chain.fasta" -x "*.debug" -x "*.valid" tmp/tl2q1 tmp/tr2${s}${i}
	
		sed -n 0~2p tmp/tr2${s}${i}/s2_chain.fasta | sort > s2sorted.txt
		diff tl2_s2_chains.sorted.txt s2sorted.txt
			sort tmp/tr2${s}${i}/s2_compressed_chain.debug > s2sortedc.txt
			diff tl2_s2_chains.sorted.txt s2sortedc.txt
	
		sed -n 0~2p tmp/tr2${s}${i}/s_chain.fasta | sort > ssorted.txt
		diff tl2_s_chains.sorted.txt ssorted.txt
			sort tmp/tr2${s}${i}/s_compressed_chain.debug > ssortedc.txt
			diff tl2_s_chains.sorted.txt ssortedc.txt

		
	done



done
