# POP
The final project in Bioinformatics Institute.

Student: Dmitrii Podgalo

Supervisor: Petr Popov, Skoltech.

## Structure-based modeling of cysteine and serine disease variants of human proteome

### Motivation
There are many disease-associated mutations that endow pharmacological target (typically a protein), drug resistance, e.g. G12C amino acid substitution in oncogenic target KRAS.
People carrying such mutations may need the development of personalized drugs, that take into account structural peculiarities of the mutated protein.
One of the promising strategies is to develop covalent drugs that are specific for a given mutation.

### Goal
The goal of this project is to model structures of human proteins with disease-associated amino acid substitutions. Two types of amino acid substitutions are selected: X to Cysteine or X to Serine (X is any amino acid residue) &ndash; these residues are often used as the attachment points for covalent drugs.

### Outcome
The obtained structural models will be used as the starting conformations for the structure-based drug design pipelines.

## Description
It is require only python 3 and conda. All the work is done with command line using **bash**. This project use **GRCh37/hg19** human assembly.

There are 5 sections:

- [Statistics](#STAT) - some facts from results,
- [Installation](#INSTALL) - what to install before work, 
- [Downloads](#DOWNLOAD) - what to download before work,
- [Laboratory journal (pipeline)](#LAB) - step-by-step description of this work,
- [References](#REF).

Normal models are in **normal_models** folder. 

Mutant models are in **mutant_models** folder.

**sum.tsv** contains the project summary:

1. protein name
2. uniprot id
3. normal sequence
4. mutant position
5. normal amino acid
6. mutant amino acid
7. AlphaFold2 score for this position
8. associated disease
9. link to AlphaFold2 model (look at **normal_models** folder)

**sum.tsv head:**
```
HFM1_HUMAN	A2PYH4  <SEQUENCE>  884 I	S	94.26	Premature_ovarian_failure_9 AF-A2PYH4-F1-model_v2.pdb.gz
TM218_HUMAN	A2RU14	<SEQUENCE>	80	R	C	88.71	Joubert_syndrome	AF-A2RU14-F1-model_v2.pdb.gz
S38A8_HUMAN	A6NNN8	<SEQUENCE>	32	I	S	89.86	Foveal_hypoplasia	AF-A6NNN8-F1-model_v2.pdb.gz
GRCR1_HUMAN	A8MXD5	<SEQUENCE>	138	R	C	90.27	Autosomal_recessive_nonsyndromic_hearing_loss_25	AF-A8MXD5-F1-model_v2.pdb.gz
KBTBD_HUMAN	C9JR72	<SEQUENCE>	248	R	S	97.59	Nemaline_myopathy_6	AF-C9JR72-F1-model_v2.pdb.gz
```

<a name="INSTALL"/>

## Statistics
Some statistics from result file **sum.tsv**.

<img src="https://raw.github.com/DmitriiPodgalo/POP/main/1.png" alt="drawing" width="500"/>

Proteins distribution

<img src="https://raw.github.com/DmitriiPodgalo/POP/main/2.png" alt="drawing" width="500"/>

Normal amino acid distribution

<img src="https://raw.github.com/DmitriiPodgalo/POP/main/3.png" alt="drawing" width="500"/>

Mutant amino acid distribution

<img src="https://raw.github.com/DmitriiPodgalo/POP/main/4.png" alt="drawing" width="500"/>

Diseases distribution


<a name="INSTALL"/>

## Installation

**Install VEP for command line:**
```
conda install -c bioconda ensembl-vep
```
**Download Rosetta from (lisence required):**
https://www.rosettacommons.org/software/license-and-download

**For Rosetta need SCons**:
https://scons.org

**Get to scons folder and run:**
```
python setup.py install
```

**Then tar rosetta archive:**
```
tar -xvzf rosetta.tgz
```

**Get to rosetta path and run:**
```
path/to/rosetta/main/source/scons.py -j40 mode=release extras=cxx11thread bin
```

<a name="DOWNLOAD"/>

## Downloads

**Get gnomAD 2.1.1 ExAC exomes:**
```
wget -c -O exac.vcf.gz https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/ExAC.r1.sites.vep.vcf.gz
```
**exac.vcf.gz head:**
```
1       13372   .       G       C       608.91  PASS    <COLUMN WITH ADDITIONAL INFO>
1       13380   .       C       G       7829.15 VQSRTrancheSNP99.60to99.80  <COLUMN WITH ADDITIONAL INFO>
1       13382   .       C       G       320.40  VQSRTrancheSNP99.60to99.80  <COLUMN WITH ADDITIONAL INFO>
1       13402   .       G       C       89.66   VQSRTrancheSNP99.60to99.80  <COLUMN WITH ADDITIONAL INFO>
1       13417   .       C       CGAGA   258189.04       PASS    <COLUMN WITH ADDITIONAL INFO>
```

**Get human database for VEP:**
```
cd $HOME/.vep
wget -с http://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_vep_105_GRCh37.tar.gz
tar xzf homo_sapiens_vep_105_GRCh37.tar.gz
```

**Get clinvar database for VEP:**
```
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
```
**clinvar.vcf.gz head:**
```
1       861332  1019397 G       A       .       .   <DISEASE INFO>
1       861336  1543320 C       T       .       .   <DISEASE INFO>
1       861349  1648427 C       T       .       .   <DISEASE INFO>
1       861356  1362713 T       C       .       .   <DISEASE INFO>
1       861366  1568423 C       T       .       .   <DISEASE INFO>
```

**Get AlphaFold2 database:**
```
wget -c https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v2.tar
```

<a name="LAB"/>

## Laboratory journal (pipeline)

### Run VEP
**Run VEP to annotate:**
```
vep --input_file exac.vcf.gz \
    --output_file exac.vep.gz \
    --compress_output gzip \ # compress results
    --offline \              # use local database
    --GRCh37 \               # use GRCh37 human assembly
    --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNSIGCONF,CLNSIGINCL,CLNVC,RS \ # use clinvar
    --uniprot \              # add uniprot id
    --fork 40                # parallel execution
```
**exac.vep.gz head:**
```
1_13372_G/C     1:13372 C       ENSG00000223972 ENST00000456328 Transcript      non_coding_transcript_exon_variant      620     -       -       -       -       -       IMPACT=MODIFIER;STRAND=1
1_13372_G/C     1:13372 C       ENSG00000227232 ENST00000488147 Transcript      downstream_gene_variant -       -       -       -       -       -       IMPACT=MODIFIER;DISTANCE=1032;STRAND=-1
1_13372_G/C     1:13372 C       ENSG00000223972 ENST00000515242 Transcript      non_coding_transcript_exon_variant      613     -       -       -       -       -       IMPACT=MODIFIER;STRAND=1
1_13372_G/C     1:13372 C       ENSG00000223972 ENST00000518655 Transcript      intron_variant,non_coding_transcript_variant    -       -       -       -       -       -       IMPACT=MODIFIER;STRAND=1
1_13372_G/C     1:13372 C       ENSG00000227232 ENST00000538476 Transcript      downstream_gene_variant -       -       -       -       -       -       IMPACT=MODIFIER;DISTANCE=1039;STRAND=-1
```

**Run VEP for filter:**
```
filter_vep --input_file exac.vep.gz \
           --output_file exac.filt.vep \
           --filter "Consequence is missense_variant and (Amino_acids matches /S or Amino_acids matches /C)"
```

**exac.filt.vep head:**
```
1_69241_C/T     1:69241 T       ENSG00000186092 ENST00000335137 Transcript      missense_variant        151     151     51      P/S     Ccc/Tcc -       IMPACT=MODERATE;STRAND=1;SWISSPROT=OR4F5_HUMAN;UNIPARC=UPI0000041BC1
rs140739101     1:69428 G       ENSG00000186092 ENST00000335137 Transcript      missense_variant        338     338     113     F/C     tTt/tGt -       IMPACT=MODERATE;STRAND=1;SWISSPROT=OR4F5_HUMAN;UNIPARC=UPI0000041BC1
rs150690004     1:69496 A       ENSG00000186092 ENST00000335137 Transcript      missense_variant        406     406     136     G/S     Ggc/Agc -       IMPACT=MODERATE;STRAND=1;SWISSPROT=OR4F5_HUMAN;UNIPARC=UPI0000041BC1
1_69523_G/T     1:69523 T       ENSG00000186092 ENST00000335137 Transcript      missense_variant        433     433     145     G/C     Ggc/Tgc -       IMPACT=MODERATE;STRAND=1;SWISSPROT=OR4F5_HUMAN;UNIPARC=UPI0000041BC1
1_69634_T/A     1:69634 A       ENSG00000186092 ENST00000335137 Transcript      missense_variant        544     544     182     C/S     Tgt/Agt -       IMPACT=MODERATE;STRAND=1;SWISSPROT=OR4F5_HUMAN;UNIPARC=UPI0000041BC1
```

We chose only missense variant with substitution to serine or cysteine amino acid.

### Hand filter
**Choose variants that 100% associated with disease:**
```
cat exac.filt.vep | grep -v '^##' | grep 'CLNSIG=Pathogenic' | cut --complement -f 6,7,13 > exac.filt.pat.columns.vep
```

**Choose amino acid columns and their poristion in gene:**
```
cat exac.filt.pat.columns.vep | cut -f 8,9 > temp_1.tsv
cat temp_1.tsv | cut -f 2 | awk '{FS="/";OFS="\t"} {print $1,$2}' > AC.tsv
cat temp_1.tsv | cut -f 1 > pos.tsv
rm temp_1.tsv
```

**Choose disease and uniprac id:**
```
cat exac.filt.pat.columns.vep | cut -f 11 | awk '{FS=";";OFS="\t"} {print $0}' | grep -o 'ClinVar_CLNDN=\w*' | awk '{FS="="}{print $2}' > diseases.tsv
cat exac.filt.pat.columns.vep | cut -f 11 | awk '{FS=";";OFS="\t"} {print $0}' | grep -o 'UNIPARC=\w*' | awk '{FS="="}{print $2}' > uniprac.tsv
```

**Merge files:**
```
paste uniprac.tsv pos.tsv AC.tsv score.tsv diseases.tsv > uniprac_pos_AC_score_diseases.tsv
```

**uniprac_pos_AC_diseases.tsv head:**
```
UPI00000000ED   334     R       C    96.17      Neuronal_ceroid_lipofuscinosis_3
UPI00000000ED   334     R       C    96.17      Neuronal_ceroid_lipofuscinosis_3
UPI00000000ED   334     R       C    96.17      Neuronal_ceroid_lipofuscinosis_3
UPI000000014B   264     P       S    96.14      Blue_color_blindness
UPI0000000239   120     P       S    98.71      alpha_Thalassemia
```

### Mapping
Uniprac ID are different from AlphaFold2 ID, so we need map these ID at https://www.uniprot.org/uploadlists/. 

Just upload uniprac.tsv and choose next features: 
- UniProt, 
- ProteinName, 
- UniParc, 
- Sequence. 

Then download to uniprot.tsv.

### Merge data
**Get AlphaFold2 models list of id and path:**
```
ls path/to/alpha_fold/ | grep -v 'cif.gz' | grep -v 'UP000005640' | awk '{FS="-"}{print $2}' > id.tsv
ls ../alpha_fold/ | grep -v 'cif.gz' | grep -v 'UP000005640' > path.tsv
paste id.tsv path.tsv > id_path.tsv
```

**id_path.tsv head:**
```
A0A024R1R8      AF-A0A024R1R8-F1-model_v2.pdb.gz
A0A024RBG1      AF-A0A024RBG1-F1-model_v2.pdb.gz
A0A024RCN7      AF-A0A024RCN7-F1-model_v2.pdb.gz
A0A075B6H5      AF-A0A075B6H5-F1-model_v2.pdb.gz
A0A075B6H7      AF-A0A075B6H7-F1-model_v2.pdb.gz
```

**Inner join uniprot id and AlphaFold2 id and path:**
```
sort -k1,1 -o uniprot.tsv uniprot.tsv
sort -k1,1 -o id_path.tsv id_path.tsv
join id_path.tsv uniprot.tsv > uniprot_id_path.tsv
```

**uniprot_id_path.tsv head:**
```
P16671  AF-P16671-F1-model_v2.pdb.gz    CD36_HUMAN      UPI0000000C91   <SEQUENCE>
P51659  AF-P51659-F1-model_v2.pdb.gz    DHB4_HUMAN      UPI0000000C4F   <SEQUENCE>
Q15582  AF-Q15582-F1-model_v2.pdb.gz    BGH3_HUMAN      UPI0000000C6A   <SEQUENCE>
P15559  AF-P15559-F1-model_v2.pdb.gz    NQO1_HUMAN      UPI0000000C86   <SEQUENCE>
P16671  AF-P16671-F1-model_v2.pdb.gz    CD36_HUMAN      UPI0000000C91   <SEQUENCE>
```

**Choose columns and set delimiter as "\t":**
```
cat uniprot_id_path.tsv | awk '{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5}' > temp.tsv && mv temp.tsv uniprot_id_path.tsv
```

**Inner join of files to summary file:**
```
sort -k1,1 -o uniprac_pos_AC_diseases.tsv uniprac_pos_AC_score_diseases.tsv
sort -k1,1 -o uniprot_id_path.tsv uniprot_id_path.tsv
join -1 3 -2 1 uniprot_id_path.tsv uniprac_pos_AC_score_diseases.tsv | awk '{print $3, $2, $4, $5, $6, $7, $8, $9}' | uniq > sum.tsv
```

**Remove variants that have not disease:**
```
cat sum.tsv | grep -v 'not_\w*' > temp.tsv | mv temp.tsv sum.tsv
```

### Rosetta ddg_monomer
Used .sh script:
```
#!/bin/bash

pre_mut='total 1
1
'

for i in $(seq 1 1339);
do
        echo "$pre_mut"`cat sum.tsv | awk '{print $5, $4, $6}' | sed "$i!d"` > mut_file.txt # write mut_file
        path='path/to/alpha_fold/'`cut -f 8 sum.tsv | sed "$i!d"`                           # prepare path to normal model
        name=`cut -f 2 sum.tsv | sed "$i!d"`'_mut.pdb'                                      # name for saving
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
```

### Summary
Script to collect normal models to **normal_models** folder:
```
#!/bin/bash

for i in $(seq 1 1339);
do
        path='path/to/alpha_fold/'`cut -f 8 sum.tsv | sed "$i!d"` # prepare path to normal model
        echo $path                                                # logging
        `cp $path ./norm_models/`                                 # copy the model to the folder
done
```

Script to collect AlphaFold2 score to score.tsv:
```
#!/bin/bash

for i in $(seq 1 1339);
do
        path='../../../alpha_fold/'`cut -f 8 ../sum.tsv | sed "$i!d"`
        pos=`cat ../sum.tsv | cut -f 4 | sed "$i!d"`
        echo $path
        echo `zcat $path | grep '^ATOM' | awk '{print $6, $11}' | uniq | grep "^$pos " | cut -f 2` >> score.tsv
done
```

### Conclusion
In the begging we had 9362318 variants, after filter we got 1339 variants 100% associated with disease. For each of 1339 variants got mutant 3D-model. 

<a name="REF"/>

## References
1. Barlow, Kyle A., et al. "Flex ddG: Rosetta ensemble-based estimation of changes in protein–protein binding affinity upon mutation." The Journal of Physical Chemistry B 122.21 (2018): 5389-5399.
2. Jumper, John, et al. "AlphaFold 2." (2020).
3. Knight, Steven. "Building software with SCons." Computing in Science & Engineering 7.1 (2005): 79-88.
4. Koch, Linda. "Exploring human genomic diversity with gnomAD." Nature Reviews Genetics 21.8 (2020): 448-448.
5. Landrum, Melissa J., et al. "ClinVar: public archive of interpretations of clinically relevant variants." Nucleic acids research 44.D1 (2016): D862-D868.
6. McCarthy, Davis J., et al. "Choice of transcripts and software has a large effect on variant annotation." Genome medicine 6.3 (2014): 1-16.
