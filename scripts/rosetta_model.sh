#!/bin/bash

pre_mut='total 1
1
'

for i in $(seq 1 1339);
do
	echo "$pre_mut"`cat sum.tsv | awk '{print $5, $4, $6}' | sed "$i!d"` > mut_file.txt                                # write mut_file
	path='path/to/alpha_fold/'`cut -f 8 sum.tsv | sed "$i!d"`                                                          # prepare path to normal model
	name=`cut -f 2 sum.tsv | sed "$i!d"`'_mut_'`cat sum.tsv | awk '{OFS=""}{print $5, $4, $6}' | sed "$i!d"`'.pdb'     # name for saving
	echo $path
	`../rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/ddg_monomer.cxx11thread.linuxgccrelease -in:file:s $path -ddg::mut_file mut_file.txt -multithreading:total_threads 40` # execute rosetta
	`rm repacked*`
	`rm *1.pdb *2.pdb`
	`rm *out mutant_traj* wt_traj`
	`mv *3.pdb mut_models/$name`
	`gzip mut_models/$name`
	echo "Model $i done"
		# above - remove draft files and move model to the folder
done
