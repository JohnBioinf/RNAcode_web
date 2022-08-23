#!/bin/bash

set -e

shellcheck disable=SC1091
source ./auxilary_scripts/activate_env_backend.sh
python ./system_setup/update_ete3_ncbi_db.py

echo "Update blast database"

blast_db="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['blast_db'])" \
	< "./parameters_backend_local.json")"

blast_bin="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['blast_bin'])" \
	< "./parameters_backend_local.json")"

PATH="$blast_bin:$PATH"

if [ -z "$blast_db" ]; then
	echo "Parameter blast_db not set please edit parameter_backend.json"
	exit 1
fi

if [ ! -d "$blast_db" ]; then
	mkdir "$blast_db"
fi

# shellcheck disable=SC2164
cd "$blast_db"

update_blastdb.pl --decompress ref_euk_rep_genomes --quiet
update_blastdb.pl --decompress ref_prok_rep_genomes --quiet
update_blastdb.pl --decompress ref_viroids_rep_genomes --quiet
update_blastdb.pl --decompress ref_viruses_rep_genomes --quiet
update_blastdb.pl --decompress nt --quiet

blastdbcmd -db ref_euk_rep_genomes -entry all -outfmt %T > ref_euk_rep_genomes.taxidlist
blastdbcmd -db ref_prok_rep_genomes -entry all -outfmt %T > ref_prok_rep_genomes.taxidlist
blastdbcmd -db ref_viroids_rep_genomes -entry all -outfmt %T > ref_viroids_rep_genomes.taxidlist
blastdbcmd -db ref_viruses_rep_genomes -entry all -outfmt %T > ref_viruses_rep_genomes.taxidlist
blastdbcmd -db nt -entry all -outfmt %T > nt.taxidlist

echo "Finished"
