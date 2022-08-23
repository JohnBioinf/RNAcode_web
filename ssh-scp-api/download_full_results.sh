#!/bin/bash

set -e

user="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user'])" \
	< "./parameters_backend_local.json")"
machine_name="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['machine_name'])" \
	< "./parameters_backend_local.json")"

job_id=$1
work_dir=$2
result_dir=$3

# shellcheck disable=SC2029
ssh "$user"@"$machine_name" "tar zcvf - -C $work_dir $job_id 2> /dev/null |\
       	cat" > "$result_dir/$job_id"_full_results.tgz
