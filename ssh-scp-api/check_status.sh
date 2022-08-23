#!/bin/bash

set -e

job_status_file="$1"

user="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user'])" \
	< "./parameters_backend_local.json")"
machine_name="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['machine_name'])" \
	< "./parameters_backend_local.json")"

# shellcheck disable=SC2029
ssh "$user"@"$machine_name" cat "$job_status_file"
