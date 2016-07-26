#!/bin/sh

prefix=$1

# get the first 31mer

grep -v ">" $1_chain.fasta | cut -b-31 > $1_fasta.first31

cut -b-31 $1_chain.edges > $1_edges.first31

cut -b33-63 $1_chain.components | sort -u -T . > $1_comp.rep31

#compare
echo "compare first 31mer"

diff -w $1_fasta.first31 $1_edges.first31
diff -w $1_edges.first31 $1_comp.rep31


# get the last 31mer

grep -v ">" $1_chain.fasta | grep -o ".\{31\}$" > $1_fasta.last31

cut -b33-63 $1_chain.edges > $1_edges.last31


# compare 
echo "compare last 31mer"

diff -w $1_fasta.last31 $1_edges.last31
