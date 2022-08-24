#!/bin/bash

set -e

production_dir="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['production_dir'])" \
	< "./parameters_frontend_local.json")"
if [ -z "$production_dir" ]; then
	echo "set variable"
	exit 1
fi

yesexpr="Y|y|^$"
noexpr="N|n"
while true; do
	read -rp "Delete logs in dir:$production_dir/job_submissions (Y/n)? " yn
	if [[ "$yn" =~ $yesexpr ]]; then break; fi
	if [[ "$yn" =~ $noexpr ]]; then exit; fi
	echo "Answer Y/n."
done

rm -rf "$production_dir/job_submissions"
mkdir "$production_dir/job_submissions"
