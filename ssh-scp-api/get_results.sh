#!/bin/bash

set -e

user="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user'])" \
	< "./parameters_backend_local.json")"

machine_name="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['machine_name'])" \
	< "./parameters_backend_local.json")"

user_bin_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user_bin_path'])" \
	< "./parameters_frontend_local.json")"

PATH="$user_bin_path:$PATH"

job_id=$1
RNAcode_result_path=$2
eps_path=$3
align_fasta_path=$4
align_plot_path=$5
result_dir=$6

if ssh -q "$user"@"$machine_name" [[ -f "$RNAcode_result_path" ]]; then
	scp "$user"@"$machine_name":"$RNAcode_result_path" "$result_dir/$job_id"_rnacode_result.tsv
fi

scp "$user"@"$machine_name":"$align_fasta_path" "$result_dir/$job_id"_full.aln

# some awk to make fasta
awk -F '[[:space:]][[:space:]]+' \
	'BEGIN{i=0}\
	{if(NR!=1){if(!($1 in h)){a[i++]=$1};h[$1] = h[$1] $2;}}\
	END{for(i in a){if(a[i] == ""){continue};print(">" a[i] "\n" h[a[i]])}}' \
	"$result_dir/$job_id"_full.aln > "$result_dir/$job_id"_aln.fasta

# make alignment html
timeout 120 mview -in fasta -html data -bold -coloring \
	CLUSTAL "$result_dir/$job_id"_aln.fasta -conservation on > \
	"$result_dir/$job_id"_aln.html

if ssh -q "$user"@"$machine_name" [[ -f "$align_plot_path" ]]; then
	scp "$user"@"$machine_name":"$align_plot_path" "$result_dir/$job_id"_align.svg
fi

# shellcheck disable=SC2029
readarray -t eps_files < <(ssh "$user@$machine_name" ls "$eps_path")
for eps_file in "${eps_files[@]}"; do 
	scp "$user@$machine_name:$eps_path/$eps_file" "$result_dir/${job_id}_$eps_file"
done
