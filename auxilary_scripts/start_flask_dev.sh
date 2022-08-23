#!/bin/bash

export FLASK_DEBUG=1

python_env_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['python_env_path'])" \
	< "./parameters_frontend_local.json")"

ip="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['ip'])" \
	< "./parameters_frontend_local.json")"

port="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['port'])" \
	< "./parameters_frontend_local.json")"

# shellcheck disable=SC1091
source "$python_env_path/bin/activate"

flask run --host="$ip" --port "$port"
