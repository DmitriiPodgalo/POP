#!/bin/bash

for i in $(seq 1 1339);
do
	        path='path/to/alpha_fold/'`cut -f 8 sum.tsv | sed "$i!d"` # prepare path to normal model
		        echo $path                                                # logging
			        `cp $path ./norm_models/`                                 # copy the model to the folder
			done
