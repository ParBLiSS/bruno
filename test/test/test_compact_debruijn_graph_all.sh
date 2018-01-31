#!/bin/sh

DIR=test_dbg

rm -rf ${DIR}

#gold
mkdir -p ${DIR}/gold-q1
echo "gold q1"
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -O ${DIR}/gold-q1/s ~/src/bruno/test/data/test.star.unitiqs.fastq >   ${DIR}/gold-q1_s.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -O ${DIR}/gold-q1/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/gold-q1_s2.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -O ${DIR}/gold-q1/u ~/src/bruno/test/data/test.unitiqs.fastq > 	   ${DIR}/gold-q1_u.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -O ${DIR}/gold-q1/t ~/src/bruno/test/data/test.short.unitiqs.fastq >  ${DIR}/gold-q1_t.log 2>&1

mkdir -p ${DIR}/gold-q1-tl2
echo "gold q1 tl2"
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -T -L 2 -O ${DIR}/gold-q1-tl2/s ~/src/bruno/test/data/test.star.unitiqs.fastq > 	${DIR}/gold-q1-tl2_s.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -T -L 2 -O ${DIR}/gold-q1-tl2/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > 	${DIR}/gold-q1-tl2_s2.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -T -L 2 -O ${DIR}/gold-q1-tl2/u ~/src/bruno/test/data/test.unitiqs.fastq > 		${DIR}/gold-q1-tl2_u.log 2>&1
/usr/bin/time -v bin/compact_debruijn_graph_fastq_refactored -T -L 2 -O ${DIR}/gold-q1-tl2/t ~/src/bruno/test/data/test.short.unitiqs.fastq > 	${DIR}/gold-q1-tl2_t.log 2>&1


for i in 1 2 4 8 16
do
# no dna5
	for s in 4 16
	do 
	# blocked
		mkdir -p ${DIR}/b${s}-${i}
		echo "$i $s b"
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -C -O ${DIR}/b${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/b${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -C -O ${DIR}/b${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/b${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -C -O ${DIR}/b${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/b${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -C -O ${DIR}/b${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/b${s}-${i}_t.log 2>&1
		
	# freq

		mkdir -p ${DIR}/f${s}-${i}
		echo "$i $s f"
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -C -O ${DIR}/f${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/f${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -C -O ${DIR}/f${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/f${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -C -O ${DIR}/f${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/f${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -C -O ${DIR}/f${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/f${s}-${i}_t.log 2>&1
		
	#blocked min mem

		mkdir -p ${DIR}/bm${s}-${i}
		echo "$i $s bm"
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -C -O ${DIR}/bm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/bm${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -C -O ${DIR}/bm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/bm${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -C -O ${DIR}/bm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/bm${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -C -O ${DIR}/bm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/bm${s}-${i}_t.log 2>&1
		
	# freq min mem

		mkdir -p ${DIR}/fm${s}-${i}
		echo "$i $s fm"
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -C -O ${DIR}/fm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/fm${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -C -O ${DIR}/fm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/fm${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -C -O ${DIR}/fm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/fm${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -C -O ${DIR}/fm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/fm${s}-${i}_t.log 2>&1
		

	#blocked min mem, no opt

		mkdir -p ${DIR}/obm${s}-${i}
		echo "$i $s obm"
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -C -O ${DIR}/obm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/obm${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -C -O ${DIR}/obm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/obm${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -C -O ${DIR}/obm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/obm${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -C -O ${DIR}/obm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/obm${s}-${i}_t.log 2>&1
		
	# freq min mem, no opt

		mkdir -p ${DIR}/ofm${s}-${i}
		echo "$i $s ofm"
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -C -O ${DIR}/ofm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/ofm${s}-${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -C -O ${DIR}/ofm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/ofm${s}-${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -C -O ${DIR}/ofm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/ofm${s}-${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -C -O ${DIR}/ofm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/ofm${s}-${i}_t.log 2>&1

	done
done

for t in 1 2
do

	for i in 1 2 4 8 16
	do
	# no dna5
		for s in 4 16
		do 
			echo "$i $s"
		# blocked
			mkdir -p ${DIR}/tl${t}-b${s}-${i}
			echo "$i $s b tl${t}"
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L ${t} -C -O ${DIR}/tl${t}-b${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/tl${t}-b${s}-${i}_s.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L ${t} -C -O ${DIR}/tl${t}-b${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/tl${t}-b${s}-${i}_s2.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L ${t} -C -O ${DIR}/tl${t}-b${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/tl${t}-b${s}-${i}_u.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31 -R -T -L ${t} -C -O ${DIR}/tl${t}-b${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/tl${t}-b${s}-${i}_t.log 2>&1

		# freq
			
			mkdir -p ${DIR}/tl${t}-f${s}-${i}
			echo "$i $s f tl${t}"
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -T -L ${t} -C -O ${DIR}/tl${t}-f${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/tl${t}-f${s}-${i}_s.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -T -L ${t} -C -O ${DIR}/tl${t}-f${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/tl${t}-f${s}-${i}_s2.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -T -L ${t} -C -O ${DIR}/tl${t}-f${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/tl${t}-f${s}-${i}_u.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq -R -T -L ${t} -C -O ${DIR}/tl${t}-f${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/tl${t}-f${s}-${i}_t.log 2>&1

		#blocked min mem

			mkdir -p ${DIR}/tl${t}-bm${s}-${i}
			echo "$i $s bn tl${t}"
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-bm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/tl${t}-bm${s}-${i}_s.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-bm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/tl${t}-bm${s}-${i}_s2.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-bm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/tl${t}-bm${s}-${i}_u.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-bm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/tl${t}-bm${s}-${i}_t.log 2>&1
			
		# freq min mem
			
			mkdir -p ${DIR}/tl${t}-fm${s}-${i}
			echo "$i $s fm tl${t}"	
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-fm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/tl${t}-fm${s}-${i}_s.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-fm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/tl${t}-fm${s}-${i}_s2.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-fm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/tl${t}-fm${s}-${i}_u.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -R -T -L ${t} -C -O ${DIR}/tl${t}-fm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/tl${t}-fm${s}-${i}_t.log 2>&1

		#blocked min mem, no opt

			mkdir -p ${DIR}/tl${t}-obm${s}-${i}
			echo "$i $s obm tl${t}"
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -T -L ${t} -C -O ${DIR}/tl${t}-obm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/tl${t}-obm${s}-${i}_s.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -T -L ${t} -C -O ${DIR}/tl${t}-obm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/tl${t}-obm${s}-${i}_s2.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -T -L ${t} -C -O ${DIR}/tl${t}-obm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/tl${t}-obm${s}-${i}_u.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_incr -T -L ${t} -C -O ${DIR}/tl${t}-obm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/tl${t}-obm${s}-${i}_t.log 2>&1
		
		# freq min mem, no opt
			
			mkdir -p ${DIR}/tl${t}-ofm${s}-${i}
			echo "$i $s ofm tl${t}"
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -T -L ${t} -C -O ${DIR}/tl${t}-ofm${s}-${i}/s ~/src/bruno/test/data/test.star.unitiqs.fastq   > ${DIR}/tl${t}-ofm${s}-${i}_s.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -T -L ${t} -C -O ${DIR}/tl${t}-ofm${s}-${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fastq > ${DIR}/tl${t}-ofm${s}-${i}_s2.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -T -L ${t} -C -O ${DIR}/tl${t}-ofm${s}-${i}/u ~/src/bruno/test/data/test.unitiqs.fastq        > ${DIR}/tl${t}-ofm${s}-${i}_u.log 2>&1
			/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fastq_A${s}_K31_freq_incr -T -L ${t} -C -O ${DIR}/tl${t}-ofm${s}-${i}/t ~/src/bruno/test/data/test.short.unitiqs.fastq  > ${DIR}/tl${t}-ofm${s}-${i}_t.log 2>&1

		done
	done

done

quiet=-q

echo "compare unitiqs to original"
diff ${DIR}/gold-q1/u_chain.fasta ~/src/bruno/test/data/test.unitiqs.fasta


sed -n 0~2p ${DIR}/gold-q1/s2_chain.fasta | sort > s2_chains.sorted.txt
sed -n 0~2p ${DIR}/gold-q1/s_chain.fasta | sort > s_chains.sorted.txt
sed -n 0~2p ${DIR}/gold-q1/u_chain.fasta | sort > u_chains.sorted.txt
sed -n 0~2p ${DIR}/gold-q1/t_chain.fasta | sort > t_chains.sorted.txt

sed -n 0~2p ${DIR}/gold-q1-tl2/s2_chain.fasta | sort > tl2_s2_chains.sorted.txt
sed -n 0~2p ${DIR}/gold-q1-tl2/s_chain.fasta | sort  > tl2_s_chains.sorted.txt

echo "COMPARE TO 1 proc, refactored, filtered, fastq"

for i in 1 2 4 8 16
do
# no DNA5
	for s in 4
	do 

		for prefix in b f bm fm obm ofm
		do

			echo "compare ${prefix}${s}-${i}"
			diff ${quiet} -x "*.fasta" -x "*.debug" ${DIR}/gold-q1 ${DIR}/${prefix}${s}-${i}

			sed -n 0~2p ${DIR}/${prefix}${s}-${i}/s2_chain.fasta | sort > s2sorted.txt
			diff ${quiet} s2_chains.sorted.txt s2sorted.txt
			sort ${DIR}/${prefix}${s}-${i}/s2_compressed_chain.debug > s2sortedc.txt
			diff ${quiet} s2_chains.sorted.txt s2sortedc.txt

			sed -n 0~2p ${DIR}/${prefix}${s}-${i}/s_chain.fasta | sort > ssorted.txt
			diff ${quiet} s_chains.sorted.txt ssorted.txt
			sort ${DIR}/${prefix}${s}-${i}/s_compressed_chain.debug > ssortedc.txt
			diff ${quiet} s_chains.sorted.txt ssortedc.txt

			sed -n 0~2p ${DIR}/${prefix}${s}-${i}/u_chain.fasta | sort > usorted.txt
			diff ${quiet} u_chains.sorted.txt usorted.txt
			sort ${DIR}/${prefix}${s}-${i}/u_compressed_chain.debug > usortedc.txt
			diff ${quiet} u_chains.sorted.txt usortedc.txt

			sed -n 0~2p ${DIR}/${prefix}${s}-${i}/t_chain.fasta | sort > tsorted.txt
			diff ${quiet} t_chains.sorted.txt tsorted.txt
			sort ${DIR}/${prefix}${s}-${i}/t_compressed_chain.debug > tsortedc.txt
			diff ${quiet} t_chains.sorted.txt tsortedc.txt

			echo "compare tl1-${prefix}${s}-${i}"
			diff ${quiet} -x "*.valid" -x "*.fasta" -x "*.debug" ${DIR}/gold-q1 ${DIR}/tl1-${prefix}${s}-${i}
		
			sed -n 0~2p ${DIR}/tl1-${prefix}${s}-${i}/s2_chain.fasta | sort > s2sorted.txt
			diff ${quiet} s2_chains.sorted.txt s2sorted.txt
				sort ${DIR}/tl1-${prefix}${s}-${i}/s2_compressed_chain.debug > s2sortedc.txt
				diff ${quiet} s2_chains.sorted.txt s2sortedc.txt
		
			sed -n 0~2p ${DIR}/tl1-${prefix}${s}-${i}/s_chain.fasta | sort > ssorted.txt
			diff ${quiet} s_chains.sorted.txt ssorted.txt
				sort ${DIR}/tl1-${prefix}${s}-${i}/s_compressed_chain.debug > ssortedc.txt
				diff ${quiet} s_chains.sorted.txt ssortedc.txt
		
			sed -n 0~2p ${DIR}/tl1-${prefix}${s}-${i}/u_chain.fasta | sort > usorted.txt
			diff ${quiet} u_chains.sorted.txt usorted.txt
				sort ${DIR}/tl1-${prefix}${s}-${i}/u_compressed_chain.debug > usortedc.txt
				diff ${quiet} u_chains.sorted.txt usortedc.txt
		
			sed -n 0~2p ${DIR}/tl1-${prefix}${s}-${i}/t_chain.fasta | sort > tsorted.txt
			diff ${quiet} t_chains.sorted.txt tsorted.txt
				sort ${DIR}/tl1-${prefix}${s}-${i}/t_compressed_chain.debug > tsortedc.txt
				diff ${quiet} t_chains.sorted.txt tsortedc.txt

			echo "compare tl2-${prefix}${s}-${i}"
			diff ${quiet} -x "*_chain.fasta" -x "*.debug" -x "*.valid" ${DIR}/gold-q1-tl2 ${DIR}/tl2-${prefix}${s}-${i}
		
			sed -n 0~2p ${DIR}/tl2-${prefix}${s}-${i}/s2_chain.fasta | sort > s2sorted.txt
			diff ${quiet} tl2_s2_chains.sorted.txt s2sorted.txt
				sort ${DIR}/tl2-${prefix}${s}-${i}/s2_compressed_chain.debug > s2sortedc.txt
				diff ${quiet} tl2_s2_chains.sorted.txt s2sortedc.txt
		
			sed -n 0~2p ${DIR}/tl2-${prefix}${s}-${i}/s_chain.fasta | sort > ssorted.txt
			diff ${quiet} tl2_s_chains.sorted.txt ssorted.txt
				sort ${DIR}/tl2-${prefix}${s}-${i}/s_compressed_chain.debug > ssortedc.txt
				diff ${quiet} tl2_s_chains.sorted.txt ssortedc.txt
		done		
	done

done
