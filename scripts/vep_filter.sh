#!/bin/bash

filter_vep --input_file exac.vep.gz \
	   --output_file exac.filt.vep \
	   --filter "Consequence is missense_variant and (Amino_acids matches /S or Amino_acids matches /C)"
