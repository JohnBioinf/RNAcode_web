#!/bin/bash

set -e

python_env_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['python_env_backend_path'])" \
	< "./parameters_backend_local.json")"

user="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user'])" \
	< "./parameters_backend_local.json")"

machine_name="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['machine_name'])" \
	< "./parameters_backend_local.json")"

RNAcode_web_repo_backend="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['RNAcode_web_repo_backend'])" \
	< "./parameters_backend_local.json")"

job_id=$1
min_pair_dist=$2
max_pair_dist=$3
genome_start=$4
db_type=$5
input_seq_nuc=$6

# shellcheck disable=SC2029
ssh "$user"@"$machine_name" "source $python_env_path/bin/activate &&\
    cd $RNAcode_web_repo_backend &&\
    nohup xvfb-run -d python3 $RNAcode_web_repo_backend/RNAcodeWebCore_NCBI_DB.py  \
    $job_id $min_pair_dist $max_pair_dist $genome_start $db_type $input_seq_nuc \
    &> /dev/null </dev/null " &
ssh_pid=$!
for i in $(seq 1 6); do
    if ssh -q "$user"@"$machine_name" [[ -f "$RNAcode_web_repo_backend/logs/$job_id.out" ]]; then
        break
    fi
    sleep 4
    if [[ $i -gt 5 ]]; then
        kill -9 "$ssh_pid" || true
        >&2 echo "Can not start job backend!"
        >&2 echo "Commands:"
        >&2 echo "ssh $user@$machine_name"
        >&2 echo "source $python_env_path/bin/activate"
        >&2 echo "cd $RNAcode_web_repo_backend"
        >&2 echo "nohup xvfb-run -d python3 $RNAcode_web_repo_backend/RNAcodeWebCore_NCBI_DB.py $job_id $min_pair_dist $max_pair_dist $genome_start $db_type $input_seq_nuc"
        exit 1
    fi
done
kill -9 "$ssh_pid" || true
