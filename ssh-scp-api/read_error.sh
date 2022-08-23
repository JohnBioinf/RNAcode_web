#!/bin/bash

set -e

job_id=$1

user="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user'])" \
	< "./parameters_backend_local.json")"

machine_name="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['machine_name'])" \
	< "./parameters_backend_local.json")"

RNAcode_web_repo_backend="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['RNAcode_web_repo_backend'])" \
	< "./parameters_backend_local.json")"

# shellcheck disable=SC2029
ssh "$user"@"$machine_name" cat "$RNAcode_web_repo_backend/logs/$job_id.out"
# shellcheck disable=SC2029
ssh "$user"@"$machine_name" cat "$RNAcode_web_repo_backend/logs/$job_id.err"
