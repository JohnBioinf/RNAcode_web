#!/bin/bash

set -e

user="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user'])" \
	< "./parameters_backend_local.json")"

machine_name="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['machine_name'])" \
	< "./parameters_backend_local.json")"

current_workdir=$1
out_file=$2
err_file=$3

# shellcheck disable=SC2029
ssh "$user"@"$machine_name" "rm -r $current_workdir"
# shellcheck disable=SC2029
ssh "$user"@"$machine_name" "rm -r $out_file"
# shellcheck disable=SC2029
ssh "$user"@"$machine_name" "rm -r $err_file"
