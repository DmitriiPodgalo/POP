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

## Into
It is require only python 3 and conda. All the work is done with command line using **bash**. This project use **GRCh37/hg19** human assembly.

There are 4 sections:

- [STATISTICS](#STAT) - some facts from results,
- [INSTALL DEPENDENCIES](#INSTALL) - what to install before work, 
- [DOWNLOAD](#DOWNLOAD) - what to download before work,
- [LABORATORY JOURNAL (pipeline)](#LAB) - step-by-step description of this work.

Normal models are in **normal_models** folder. Mutant models are in **mutant_models** folder.

<a name="INSTALL"/>

## STATISTICS

![ScreenShot](https://{https://drive.google.com/file/d/1SKGPGu5AmteuGqtbFE3xBI2g7Ojs9VOb/view?usp=sharing})


<a name="INSTALL"/>

## INSTALL DEPENDENCIES

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

## DOWNLOAD

**Get gnomAD 2.1.1 ExAC exomes:**
```
wget https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/ExAC.r1.sites.vep.vcf.gz
gunzip -c exac.vcf ExAC.r1.sites.vep.vcf.gz
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

**Get AlphaFold2 database:**
```
wget -c https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v2.tar
```

<a name="LAB"/>

## LABORATORY JOURNAL (pipeline)

### Run VEP
**Run VEP to annotate:**
```
vep --input_file exac.vcf.gz \
    --output_file exac.vep.gz \
    --compress_output gzip \
    --offline \
    --GRCh37 \
    --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNSIGCONF,CLNSIGINCL,CLNVC,RS \
    --uniprot \
    --fork 40
```

**Run VEP for filter:**
```
filter_vep --input_file exac.evth.vep.gz \
           --output_file exac.evth.filt.vep \
           --filter "Consequence is missense_variant and (Amino_acids matches /S or Amino_acids matches /C)"
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
paste uniprac.tsv pos.tsv AC.tsv diseases.tsv > uniprac_pos_AC_diseases.tsv
```

### Mapping
Uniprac ID are different from AlphaFold2 ID, so we need map these ID at https://www.uniprot.org/uploadlists/. Just upload uniprac.tsv and choose next features: UniProt, ProteinName, UniParc, Sequence. Then download to uniprot.tsv.

### Merge data
**Get AlphaFold2 models list of id and path:**
```
ls path/to/alpha_fold/ | grep -v 'cif.gz' | grep -v 'UP000005640' | awk '{FS="-"}{print $2}' > id.tsv
ls ../alpha_fold/ | grep -v 'cif.gz' | grep -v 'UP000005640' > path.tsv
paste id.tsv path.tsv > id_path.tsv
```

**Inner join uniprot id and AlphaFold2 id and path:**
```
sort -k1,1 -o uniprot.tsv uniprot.tsv
sort -k1,1 -o id_path.tsv id_path.tsv
join uniprot.tsv id_path.tsv > uniprot_id_path.tsv
```

**Choose columns and set delimiter as "\t":**
```
cat uniprot_id_path.tsv | awk '{FS=" ";OFS="\t"}{print $1,$2,$3,$4,$5}' > temp.tsv && mv temp.tsv uniprot_id_path.tsv
```

**Inner join of files to summary file:**
```
sort -k1,1 -o uniprac_pos_AC_diseases.tsv uniprac_pos_AC_diseases.tsv
sort -k1,1 -o uniprot_id_path.tsv uniprot_id_path.tsv
join -1 3 -2 1 uniprot_id_path.tsv uniprac_pos_AC_diseases.tsv | awk '{print $3, $2, $4, $5, $6, $7, $8}' | uniq > sum.tsv
```

**Remove variants that have not disease:**
```
cat sum.tsv | grep -v 'not_\w*' > temp.tsv | mv temp.tsv sum.tsv
```

### Rosetta ddg_monomer
Used .sh script:
```
pre_mut='total 1
1
'

for i in $(seq 1 1339);
do
        echo "$pre_mut"`cat sum.tsv | awk '{print $5, $4, $6}' | sed "$i!d"` > mut_file.txt
        path='path/to/alpha_fold/'`cut -f 8 sum.tsv | sed "$i!d"`
        name=`cut -f 2 sum.tsv | sed "$i!d"`'_mut.pdb'
        echo $path
        `../rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/ddg_monomer.cxx11thread.linuxgccrelease -in:file:s $path -ddg::mut_file mut_file.txt -multithreading:total_threads 40`
        `rm repacked*`
        `rm *1.pdb *2.pdb`
        `rm *out mutant_traj* wt_traj`
        `mv *3.pdb mut_models/$name`
        `gzip mut_models/$name`
        echo "Model $i done"
done
```

### Summary
Script to collect normal models:
```
for i in $(seq 1 1339);
do
        path='path/to/alpha_fold/'`cut -f 8 sum.tsv | sed "$i!d"`
        echo $path
        `cp $path ./norm_models/`
done
```

### Conclusion
In the begging we had 9362318 variants, after filter we got 1339 variants 100% associated with disease. For each of 1339 variants got mutant 3D-model. 
