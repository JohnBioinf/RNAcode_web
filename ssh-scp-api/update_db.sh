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

# shellcheck disable=SC2029
ssh "$user"@"$machine_name" "source $python_env_path/bin/activate &&\
	cd $RNAcode_web_repo_backend &&\
	. ./system_setup/update_blast_db.sh > /dev/null" &
