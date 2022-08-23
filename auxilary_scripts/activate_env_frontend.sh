#!/bin/bash

python_env_frontend_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['python_env_path'])" \
	< "./parameters_frontend_local.json")"

# shellcheck disable=SC1090
source "$python_env_frontend_path/bin/activate"
