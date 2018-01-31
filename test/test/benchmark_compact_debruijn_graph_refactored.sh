#!/bin/sh

DIR=benchmark

rm -rf ${DIR}

for i in 1 2 4 8 # 16
do
	for s in 'a' 'q'
	do 
		
		mkdir -p ${DIR}/r${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -O ${DIR}/r${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/r${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -O ${DIR}/r${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/r${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -O ${DIR}/r${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/r${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -O ${DIR}/r${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/r${s}${i}_t.log 2>&1
		
		mkdir -p ${DIR}/tr1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 1 -O ${DIR}/tr1${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tr1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 1 -O ${DIR}/tr1${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tr1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 1 -O ${DIR}/tr1${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/tr1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 1 -O ${DIR}/tr1${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tr1${s}${i}_t.log 2>&1
	
		mkdir -p ${DIR}/tr2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 2 -O ${DIR}/tr2${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tr2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 2 -O ${DIR}/tr2${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tr2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 2 -O ${DIR}/tr2${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s} 		  > ${DIR}/tr2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -T -L 2 -O ${DIR}/tr2${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tr2${s}${i}_t.log 2>&1



		mkdir -p ${DIR}/ro${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -O ${DIR}/ro${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/ro${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -O ${DIR}/ro${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/ro${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -O ${DIR}/ro${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/ro${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -O ${DIR}/ro${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/ro${s}${i}_t.log 2>&1
		
		mkdir -p ${DIR}/tro1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tro1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tro1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/tro1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tro1${s}${i}_t.log 2>&1
	
		mkdir -p ${DIR}/tro2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tro2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tro2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s} 		  > ${DIR}/tro2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_blocked -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tro2${s}${i}_t.log 2>&1





		mkdir -p ${DIR}/fr${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -O ${DIR}/r${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/fr${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -O ${DIR}/r${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/fr${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -O ${DIR}/r${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/fr${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -O ${DIR}/r${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/fr${s}${i}_t.log 2>&1
		
		mkdir -p ${DIR}/tfr1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 1 -O ${DIR}/tr1${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tfr1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 1 -O ${DIR}/tr1${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tfr1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 1 -O ${DIR}/tr1${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/tfr1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 1 -O ${DIR}/tr1${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tfr1${s}${i}_t.log 2>&1
	
		mkdir -p ${DIR}/tfr2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 2 -O ${DIR}/tr2${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tfr2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 2 -O ${DIR}/tr2${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tfr2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 2 -O ${DIR}/tr2${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s} 		  > ${DIR}/tfr2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -T -L 2 -O ${DIR}/tr2${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tfr2${s}${i}_t.log 2>&1



		mkdir -p ${DIR}/fro${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -O ${DIR}/ro${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/rfo${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -O ${DIR}/ro${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/rfo${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -O ${DIR}/ro${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/rfo${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -O ${DIR}/ro${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/rfo${s}${i}_t.log 2>&1
		
		mkdir -p ${DIR}/tfro1${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tfro1${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tfro1${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s}        > ${DIR}/tfro1${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 1 -O ${DIR}/tro1${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tfro1${s}${i}_t.log 2>&1
	
		mkdir -p ${DIR}/tfro2${s}${i}
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/s ~/src/bruno/test/data/test.star.unitiqs.fast${s}   > ${DIR}/tfro2${s}${i}_s.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/s2 ~/src/bruno/test/data/test.star2.unitiqs.fast${s} > ${DIR}/tfro2${s}${i}_s2.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/u ~/src/bruno/test/data/test.unitiqs.fast${s} 		  > ${DIR}/tfro2${s}${i}_u.log 2>&1
		/usr/bin/time -v mpirun -np $i bin/compact_debruijn_graph_fast${s}_freq -B -R -T -L 2 -O ${DIR}/tro2${s}${i}/t ~/src/bruno/test/data/test.short.unitiqs.fast${s}  > ${DIR}/tfro2${s}${i}_t.log 2>&1

		
	done

done

