# RNAcode Web Results

This directory contains the input, output and auxiliary files generated by
RNAcode Web to generate a custom blast data base. The data base will later by
used by child process.

If the program ran successfully this is how the directory structure should look
like:

```
+-- 1_blastn.result
+-- 1_blastn.sh
+-- 1_selected_sequences.fasta
+-- 1_seqSel.sh
+-- 1_taxid_list.txt
+-- blast_db
|   +-- job_id.ndb
|   +-- job_id.nhr
|   +-- job_id.nin
|   +-- job_id.nog
|   +-- job_id.nos
|   +-- job_id.not
|   +-- job_id.nsq
|   +-- job_id.ntf
|   +-- job_id.nto
+-- blast_db.fasta
+-- blastn_1.err
+-- blastn_1.out
+-- buildDB.err
+-- buildDB.out
+-- build_db.sh
+-- input.fasta
+-- job_status.json
+-- positions_list_blastdbcmd.txt
+-- python_requirements.txt
+-- README.md
+-- selected_species.txt
+-- seqSel_1.err
+-- seqSel_1.out
+-- SeqSelection_pipeline.py
+-- SeqSelection.py
+-- TaxIDMapFile
```

For each step a bash file is generated which shows the exact call.

The dependencies for both python scripts can be found in
`python_requirements.txt`.

The two steps blastn and sequence selection might run in multiple iterations
until enough sequence had been found.

For the whole pipeline the file `job_status.json` stored the progress for each
step.

## Blastn

### Version: {{ blast_version }}
### Version (DB): {{ blast_db_version }}

In the first step the input nucleotide sequence `input.fasta` is blasted against
a NCBI data base. Either the RefSeq DB "Representative Genome Database" or the
nt DB "non redundant nucleotide Database". The location of the data base is
stored in the variable `$BLASTDB`. The file `n_taxid_list.txt` is a list of all
taxids which are included in the blast search. The list contains all sequence
that belong to a certain taxonomic rank. But excludes all previously found
sequences and all taxids which are not present in the current data base.

The associated files with this step are:

 * `input.fasta`
 * `blastn_n.err`
 * `blastn_n.out`
 * `n_blastn.result`
 * `n_blastn.sh`
 * `n_taxid_list.txt`

## Sequence selection

In this step the blast results are selected such that only one sequence per
species in chosen. Further the target sequence is expanded to include possible
surrounding exons.

The associated files with this step are:

 * `n_blastn.result`
 * `n_seqSel.sh`
 * `n_selected_sequences.fasta`
 * `seqSel_n.err`
 * `seqSel_n.out`
 * `selected_species.txt`
 * `SeqSelection_pipeline.py`
 * `SeqSelection.py`

## Makeblastdb

### Version: {{ blast_version }}

In the last step all previously selected sequences are collected and a blast DB
is constructed.

The sequences were collected by the web service and are all unique sequences in
all `n_selected_sequences.fasta.` The sequence were stored in `blast_db.fasta`.

The associated files with this step are:

 * `blast_db.fasta`
 * `blast_db/*`
 * `buildDB.out`
 * `buildDB.err`
 * `TaxIDMapFile`
