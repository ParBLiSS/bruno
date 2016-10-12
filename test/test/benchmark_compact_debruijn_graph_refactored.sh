#!/bin/sh

DIR=benchmark

rm -rf ${DIR}

for i in 1 2 4 8 16
do
	for s in 'a' 'q'
	do 
		
		mkdir -p ${DIR}/r${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -O ${DIR}/r${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > ${DIR}/r${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -O ${DIR}/r${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${DIR}/r${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -O ${DIR}/r${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s}        > ${DIR}/r${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -O ${DIR}/r${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > ${DIR}/r${s}${i}_t.log 2>&1
		
		mkdir -p ${DIR}/tr1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 1 -O ${DIR}/tr1${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tr1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 1 -O ${DIR}/tr1${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tr1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 1 -O ${DIR}/tr1${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s}        > ${DIR}/tr1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 1 -O ${DIR}/tr1${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tr1${s}${i}_t.log 2>&1
	
		mkdir -p ${DIR}/tr2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 2 -O ${DIR}/tr2${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tr2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 2 -O ${DIR}/tr2${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tr2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 2 -O ${DIR}/tr2${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s} 		  > ${DIR}/tr2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -T -L 2 -O ${DIR}/tr2${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tr2${s}${i}_t.log 2>&1



		mkdir -p ${DIR}/ro${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -O ${DIR}/ro${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > ${DIR}/ro${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -O ${DIR}/ro${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${DIR}/ro${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -O ${DIR}/ro${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s}        > ${DIR}/ro${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -O ${DIR}/ro${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > ${DIR}/ro${s}${i}_t.log 2>&1
		
		mkdir -p ${DIR}/tro1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tro1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tro1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s}        > ${DIR}/tro1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tro1${s}${i}_t.log 2>&1
	
		mkdir -p ${DIR}/tro2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/s ~/src/bliss/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tro2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/s2 ~/src/bliss/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tro2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/u ~/src/bliss/test/data/test.unitiqs.fast${s} 		  > ${DIR}/tro2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_refactored_fast${s} -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/t ~/src/bliss/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tro2${s}${i}_t.log 2>&1


		
	done

done

