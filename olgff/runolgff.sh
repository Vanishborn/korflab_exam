#!/bin/bash

gff1=$1
gff2=$2
inner=$3
outer=$4

if [ -z "$gff1" ] || [ -z "$gff2" ] || [ -z "$inner" ] || [ -z "$outer" ]
then
	echo "Four arguments expected: gff file 1, gff file2, inner loop count, outer loop count."
	exit 1
fi

for j in $(seq 1 $outer)
do
	echo "Starting outer loop iteration: $j"
	for i in $(seq 1 $inner)
	do
		echo "Starting inner loop iteration: $i"
	   ./olgff.py $gff1 $gff2 -z $i -o "overlaps${i}.tsv" >> "output${j}.txt"
	done
done

echo "--Starting txt2tsv--"
mkdir zone-time-tsv
for k in $(seq 1 $outer)
do
	echo "Converting txt output number $k to tsv"
	./txt2tsv.py -i "output${k}.txt" -o "./zone-time-tsv/${gff1}_${gff2}_1to${inner}ZONES_${k}.tsv"
done

echo "--Merging TSV--"
awk 'BEGIN {OFS="\t"} {a[$1]=a[$1] OFS $2} END{for (i in a) print i a[i]}' ./zone-time-tsv/*.tsv | sort -k1,1n > ./zone-time-tsv/merged.tsv
