#!/bin/bash

vep --input_file exac.vcf.gz \
	    --output_file exac.vep.gz \
	        --compress_output gzip \ # compress results
		    --offline \              # use local database
		        --GRCh37 \               # use GRCh37 human assembly
			    --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNSIGCONF,CLNSIGINCL,CLNVC,RS \ # use clinvar
			        --uniprot \              # add uniprot id
				    --fork 40                # parallel execution
