#!/bin/bash

for i in $(seq 1 1339);
do
	        path='../../../alpha_fold/'`cut -f 8 ../sum.tsv | sed "$i!d"`
		        pos=`cat ../sum.tsv | cut -f 4 | sed "$i!d"`
			        echo $path
				        echo `zcat $path | grep '^ATOM' | awk '{print $6, $11}' | uniq | grep "^$pos " | cut -f 2` >> score.tsv
				done
