#!/bin/bash

set -e

user="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user'])" \
	< "./parameters_backend_local.json")"

job_id=$1

if [[ $job_id =~ all|^$ ]]; then
	yesexpr="Y|y|^$"
	noexpr="N|n"
	while true; do
	    read -rp "Kill all proccess (Y/n)? " yn
	    if [[ "$yn" =~ $yesexpr ]]; then break; fi
	    if [[ "$yn" =~ $noexpr ]]; then exit; fi
	    echo "Answer Y/n."
	done
	job_reg_pat=".*"
else
	job_reg_pat=" $job_id(-child_[0-9]+|) "
	job_reg_pat_slurm="$job_id?(-child_[0-9]+|)"
fi

mapfile -t pids < <(pgrep -f "RNAcodeWebCore")
for pid in "${pids[@]}"; do
	full_command=$(ps -o command "$pid" | tail -n 1)
	prog=$(echo "$full_command" | awk '{print $1}')
	if [[ ! "$prog" =~ python|bash|/usr/bin/sh ]]; then
		continue
	fi
	if [[ ! "$full_command" =~ $job_reg_pat ]]; then
		continue
	fi
	kill "$pid"
done

if [[ $job_id =~ all|^$ ]]; then
	# Also cleans up old Xvfb instances
	pkill -f Xvfb || true
    readarray -t job_ids < <(squeue -o %20i%100u | grep "$user" | awk '{print $1}')
    for job_id in "${job_ids[@]}"; do
        scancel "$job_id"
    done
else
	# This gives 100 chars for job name currently job_id lenght is max 80
    readarray -t job_ids < <(squeue -o %20i%100j | grep -E "$job_reg_pat_slurm" | awk '{print $1}')
    for job_id in "${job_ids[@]}"; do
        scancel "$job_id"
    done
fi
