#!/bin/bash

python_env_backend_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['python_env_backend_path'])" \
	< "./parameters_backend_local.json")"

# shellcheck disable=SC1090
source "$python_env_backend_path/bin/activate"
