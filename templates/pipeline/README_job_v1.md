# RNAcode Web Results

This directory contains the input, output and auxiliary files generated by
RNAcode Web.

If the program ran successfully this is how the directory structure should look
like:

```

+-- aligment.err
+-- aligment.out
+-- align.sh
+-- blastn.err
+-- blastn.out
+-- blastn.result
+-- blastn.sh
+-- candidates_gi_list.txt
+-- candidates.fasta
+-- eps
|   +-- hss-0.eps
|   +-- hss-1.eps
|   ...
+-- full.aln
+-- input.fasta
+-- job_status.json
+-- multi.fasta
+-- RNAcode.err
+-- RNAcode.out
+-- rnacode_result.tsv
+-- RNAcode.sh
+-- seqSel.err
+-- seqSel.out
+-- seq_sel.sh
+-- taxid_list.txt
+-- align.tree

```

For each step a bash file is stored which shows the exact call for each step.
To run all computations again you can execute `run_all.sh`. Before all
variables must be set and all programs must be found in the current PATH. The
software versions for each step can be found below. For each step an `.err` and
`.out` exist which stored the output and error stream for the program.

For the whole pipeline the file `job_status.json` stored the progress for each
step.

## Blastn

### Version: {{ blast_version }}
### Version (DB): {{ blast_db_version }}

In the first step the input nucleotide sequence `input.fasta` is blasted
against the non redundant blast db (nt). The location of the data base is
stored in the variable `$BLASTDB`. The file `taxid_list.txt` is a list of all
taxids which are excluded from the blast search. The list was generated by
RNAcode web and includes all know descendants of the parameter `species_taxid`
provided by the user when the job was submitted.

The associated files with this step are:

 * `input.fasta`
 * `blastn.err`
 * `blastn.out`
 * `blastn.result`
 * `blastn.sh`
 * `taxid_list.txt`

## Sequence selection

In this step the blast results are selected such that only descendants of the
parameter `genus_taxid` are included and the sequence are from different
species. Also all sequence are checked if the not exceed the maximal sequence
similarity set by the user. Further if the user chose the option to expand the
target sequence, the program checks if the target sequences are included in the
results.

The associated files with this step are:

 * `blastn.result`
 * `candidates_gi_list.txt`
 * `candidates.fasta`
 * `seqSel.err`
 * `seqSel.out`
 * `seqSel.sh`
 * `multi.fasta`

## Alignment

### Version: {{ aligner_version }}

Next all sequences are aligned using {{ name_aligner }}. The guide tree that is
used for the alignment is printed in `align.tree`.

The associated files with this step are:

 * `multi.fasta`
 * `aligment.err`
 * `aligment.out`
 * `align.sh`
 * `full.aln`
 * `align.tree`

## RNAcode

### Version: {{ rnacode_version }}

In the last step RNAcode calculates the likelihood that the target sequence is
a coding gene. The results are written in `rnacode_result.tsv` and for good
results an `.eps` plot is drawn in the director `eps`.

The associated files with this step are:

 * `full.aln`
 * `eps (dir)`
 * `RNAcode.err`
 * `RNAcode.out`
 * `rnacode_result.tsv`
 * `RNAcode.sh`

