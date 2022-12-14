# RNAcode Web Results

This directory contains the input, output and auxiliary files generated by
RNAcode Web to find coding sequence based on a previously constructed data base.
The data base should be in the super category and the exact path can be found in
`1_blastn.sh`.

If the program ran successfully this is how the directory structure should look
like:

```
+-- 1_blastn.result
+-- 1_blastn.sh
+-- 1_selected_sequences.fasta
+-- 1_seqSel.sh
+-- alignment.aln
+-- alignment.err
+-- alignment.fasta
+-- alignment.maf
+-- alignment.out
+-- alignment.sh
+-- alignment.svg
+-- alignment.tree
+-- blastn_1.err
+-- blastn_1.out
+-- candidates.fasta
+-- eps
|   +-- hss-0.eps
|   +-- hss-1.eps
+-- input.fasta
+-- job_status.json
+-- positions_list_blastdbcmd.txt
+-- python_requirements.txt
+-- RNAcode.err
+-- RNAcode.out
+-- rnacode_result.tsv
+-- RNAcode.sh
+-- seqSel_1.err
+-- seqSel_1.out
+-- SeqSelection_pipeline.py
+-- SeqSelection.py
```

For each step a bash file is generated which shows the exact call.

The dependencies for both python scripts can be found in
`python_requirements.txt`.

For the whole pipeline the file `job_status.json` stored the progress for each
step.

## Blastn

### Version: {{ blast_version }}
### Version (DB): {{ blast_db_version }}

In the first step the input nucleotide sequence `input.fasta` is blasted against
a custom blast data base, which was previously constructed. The location of the
data base is given in the call.The associated files with this step are:

 * `input.fasta`
 * `blastn_n.err`
 * `blastn_n.out`
 * `n_blastn.result`
 * `n_blastn.sh`

## Sequence selection

In the next step the sequences in the blast result will be selected to build the
alignment. Only one sequence per species is allowed, all sequences must have a
minimal distance to each other and a maximal distance to the target. Further the
sequences of the blast results will be extended, such that it covers the
complete input sequence. This is needed as most blast results do not cover the
entire query sequence.


The associated files with this step are:

 * `n_blastn.result`
 * `n_seqSel.sh`
 * `n_selected_sequences.fasta`
 * `candidates.fasta`
 * `seqSel_n.err`
 * `seqSel_n.out`
 * `selected_species.txt`
 * `SeqSelection_pipeline.py`
 * `SeqSelection.py`

## Alignment

### Version: {{ aligner_version }}

Next all sequences are aligned using {{ name_aligner }}. The guide tree that is
used for the alignment is printed in `align.tree`.

Previously all unique sequences in all `n_selected_sequences.fasta` are
collected by the web service and stored in `alignment.fasta`.

The associated files with this step are:

 * `alignment.fasta`
 * `aligment.err`
 * `aligment.out`
 * `aligment.sh`
 * `aligment.svg`
 * `aligment.aln`
 * `aligment.tree`

## RNAcode

### Version: {{ rnacode_version }}

In the last step RNAcode calculates the likelihood that the target sequence is
a coding gene. The results are written in `rnacode_result.tsv` and for good
results an `.eps` plot is drawn in the director `eps`.

The associated files with this step are:

 * `full.aln`
 * `eps/*`
 * `RNAcode.err`
 * `RNAcode.out`
 * `rnacode_result.tsv`
 * `RNAcode.sh`
