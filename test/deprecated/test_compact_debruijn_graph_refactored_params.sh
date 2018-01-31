#!/bin/sh

DIR=param_test

rm -rf ${DIR}


mkdir -p ${DIR}/q1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O ${DIR}/q1/s ~/src/bliss/test/data/test.star.unitiqs.fastq >   ${DIR}/q1_s.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O ${DIR}/q1/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > ${DIR}/q1_s2.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O ${DIR}/q1/u ~/src/bliss/test/data/test.unitiqs.fastq > 		${DIR}/q1_u.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq -O ${DIR}/q1/t ~/src/bliss/test/data/test.short.unitiqs.fastq >  ${DIR}/q1_t.log 2>&1

mkdir -p ${DIR}/tl2q1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O ${DIR}/tl2q1/s ~/src/bliss/test/data/test.star.unitiqs.fastq > 	${DIR}/tl2q1_s.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O ${DIR}/tl2q1/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > 	${DIR}/tl2q1_s2.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O ${DIR}/tl2q1/u ~/src/bliss/test/data/test.unitiqs.fastq > 		${DIR}/tl2q1_u.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_low_mem_fastq -T -L 2 -O ${DIR}/tl2q1/t ~/src/bliss/test/data/test.short.unitiqs.fastq > 	${DIR}/tl2q1_t.log 2>&1


for i in 1 2 4 8 16
do
# no dna5
	for s in 4 16
	do 

		mkdir -p ${DIR}/ro${s}-${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -O ${DIR}/ro${s}-${i}/s ~/src/bliss/test/data/test.star.unitiqs.fastq   > ${DIR}/ro${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -O ${DIR}/ro${s}-${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > ${DIR}/ro${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -O ${DIR}/ro${s}-${i}/u ~/src/bliss/test/data/test.unitiqs.fastq        > ${DIR}/ro${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -O ${DIR}/ro${s}-${i}/t ~/src/bliss/test/data/test.short.unitiqs.fastq  > ${DIR}/ro${s}-${i}_t.log 2>&1
		
		mkdir -p ${DIR}/tro1-${s}-${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 1 -O ${DIR}/tro1-${s}-${i}/s ~/src/bliss/test/data/test.star.unitiqs.fastq   > ${DIR}/tro1-${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 1 -O ${DIR}/tro1-${s}-${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > ${DIR}/tro1-${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 1 -O ${DIR}/tro1-${s}-${i}/u ~/src/bliss/test/data/test.unitiqs.fastq        > ${DIR}/tro1-${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 1 -O ${DIR}/tro1-${s}-${i}/t ~/src/bliss/test/data/test.short.unitiqs.fastq  > ${DIR}/tro1-${s}-${i}_t.log 2>&1
	
		mkdir -p ${DIR}/tro2-${s}-${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 2 -O ${DIR}/tro2-${s}-${i}/s ~/src/bliss/test/data/test.star.unitiqs.fastq   > ${DIR}/tro2-${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 2 -O ${DIR}/tro2-${s}-${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fastq > ${DIR}/tro2-${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 2 -O ${DIR}/tro2-${s}-${i}/u ~/src/bliss/test/data/test.unitiqs.fastq 		> ${DIR}/tro2-${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L 2 -O ${DIR}/tro2-${s}-${i}/t ~/src/bliss/test/data/test.short.unitiqs.fastq  > ${DIR}/tro2-${s}-${i}_t.log 2>&1
		
	done

done

quiet=-q

echo "compare unitiqs to original"
diff ${DIR}/q1/u_chain.fasta ~/src/bliss/test/data/test.unitiqs.fasta


sed -n 0~2p ${DIR}/q1/s2_chain.fasta | sort > s2_chains.sorted.txt
sed -n 0~2p ${DIR}/q1/s_chain.fasta | sort > s_chains.sorted.txt
sed -n 0~2p ${DIR}/q1/u_chain.fasta | sort > u_chains.sorted.txt
sed -n 0~2p ${DIR}/q1/t_chain.fasta | sort > t_chains.sorted.txt

sed -n 0~2p ${DIR}/tl2q1/s2_chain.fasta | sort > tl2_s2_chains.sorted.txt
sed -n 0~2p ${DIR}/tl2q1/s_chain.fasta | sort  > tl2_s_chains.sorted.txt

echo "COMPARE TO 1 proc, refactored, filtered, fastq"

for i in 1 2 4 8 16
do
# no DNA5
	for s in 4 16
	do 

		echo "compare ro${s}-${i}"
		diff ${quiet} -x "*.fasta" -x "*.debug" ${DIR}/q1 ${DIR}/ro${s}-${i}

		sed -n 0~2p ${DIR}/ro${s}-${i}/s2_chain.fasta | sort > s2sorted.txt
		diff ${quiet} s2_chains.sorted.txt s2sorted.txt
		sort ${DIR}/ro${s}-${i}/s2_compressed_chain.debug > s2sortedc.txt
		diff ${quiet} s2_chains.sorted.txt s2sortedc.txt

		sed -n 0~2p ${DIR}/ro${s}-${i}/s_chain.fasta | sort > ssorted.txt
		diff ${quiet} s_chains.sorted.txt ssorted.txt
		sort ${DIR}/ro${s}-${i}/s_compressed_chain.debug > ssortedc.txt
		diff ${quiet} s_chains.sorted.txt ssortedc.txt

		sed -n 0~2p ${DIR}/ro${s}-${i}/u_chain.fasta | sort > usorted.txt
		diff ${quiet} u_chains.sorted.txt usorted.txt
		sort ${DIR}/ro${s}-${i}/u_compressed_chain.debug > usortedc.txt
		diff ${quiet} u_chains.sorted.txt usortedc.txt

		sed -n 0~2p ${DIR}/ro${s}-${i}/t_chain.fasta | sort > tsorted.txt
		diff ${quiet} t_chains.sorted.txt tsorted.txt
		sort ${DIR}/ro${s}-${i}/t_compressed_chain.debug > tsortedc.txt
		diff ${quiet} t_chains.sorted.txt tsortedc.txt

		echo "compare tro1-${s}-${i}"
		diff ${quiet} -x "*.valid" -x "*.fasta" -x "*.debug" ${DIR}/q1 ${DIR}/tro1-${s}-${i}
	
		sed -n 0~2p ${DIR}/tro1-${s}-${i}/s2_chain.fasta | sort > s2sorted.txt
		diff ${quiet} s2_chains.sorted.txt s2sorted.txt
			sort ${DIR}/tro1-${s}-${i}/s2_compressed_chain.debug > s2sortedc.txt
			diff ${quiet} s2_chains.sorted.txt s2sortedc.txt
	
		sed -n 0~2p ${DIR}/tro1-${s}-${i}/s_chain.fasta | sort > ssorted.txt
		diff ${quiet} s_chains.sorted.txt ssorted.txt
			sort ${DIR}/tro1-${s}-${i}/s_compressed_chain.debug > ssortedc.txt
			diff ${quiet} s_chains.sorted.txt ssortedc.txt
	
		sed -n 0~2p ${DIR}/tro1-${s}-${i}/u_chain.fasta | sort > usorted.txt
		diff ${quiet} u_chains.sorted.txt usorted.txt
			sort ${DIR}/tro1-${s}-${i}/u_compressed_chain.debug > usortedc.txt
			diff ${quiet} u_chains.sorted.txt usortedc.txt
	
		sed -n 0~2p ${DIR}/tro1-${s}-${i}/t_chain.fasta | sort > tsorted.txt
		diff ${quiet} t_chains.sorted.txt tsorted.txt
			sort ${DIR}/tro1-${s}-${i}/t_compressed_chain.debug > tsortedc.txt
			diff ${quiet} t_chains.sorted.txt tsortedc.txt

		echo "compare tro2-${s}-${i}"
		diff ${quiet} -x "*_chain.fasta" -x "*.debug" -x "*.valid" ${DIR}/tl2q1 ${DIR}/tro2-${s}-${i}
	
		sed -n 0~2p ${DIR}/tro2-${s}-${i}/s2_chain.fasta | sort > s2sorted.txt
		diff ${quiet} tl2_s2_chains.sorted.txt s2sorted.txt
			sort ${DIR}/tro2-${s}-${i}/s2_compressed_chain.debug > s2sortedc.txt
			diff ${quiet} tl2_s2_chains.sorted.txt s2sortedc.txt
	
		sed -n 0~2p ${DIR}/tro2-${s}-${i}/s_chain.fasta | sort > ssorted.txt
		diff ${quiet} tl2_s_chains.sorted.txt ssorted.txt
			sort ${DIR}/tro2-${s}-${i}/s_compressed_chain.debug > ssortedc.txt
			diff ${quiet} tl2_s_chains.sorted.txt ssortedc.txt
		
	done

done
