#!/bin/bash

set -e

# Set the following variables
export BLASTDB=""
# The python enviroment can be genereated with
# python3 -m pip install -r requirements.txt
export PYTHON_ENV=""
path_blast_bin=""
path_clustal_bin=""
path_rnacode_bin=""
# set end

echo "None" > "./reference_species.txt"

path_seqSel_script="$(pwd)"
export PATH="$path_seqSel_script:$path_blast_bin:$path_rnacode_bin:$path_clustal_bin":$PATH

if [ -f ./multi.fasta ]; then
	mv ./multi.fasta ./multi.fasta_dup
fi

{%- for j in range(1, i) %}

# {{ j }}. Iteration

# call blast
bash {{ j }}_blastn.sh

# select sequences

bash {{ j }}_seqSel.sh
{%- endfor %}

# build alignment

bash alignment.sh

# RNAcode

bash RNAcode.sh

