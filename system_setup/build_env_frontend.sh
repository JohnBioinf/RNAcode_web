#!/bin/bash

set -e

R_bin="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['R_bin'])" \
	< "./parameters_frontend_local.json")"
if [ -z "$R_bin" ]; then
	echo "set variables"
	exit 1
fi

python_env_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['python_env_path'])" \
	< "./parameters_frontend_local.json")"
if [ -z "$python_env_path" ]; then
	echo "set variables"
	exit 1
fi

user_bin_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['user_bin_path'])" \
	< "./parameters_frontend_local.json")"
if [ -z "$user_bin_path" ]; then
	echo "set variables"
	exit 1
fi

production_dir="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['production_dir'])" \
	< "./parameters_frontend_local.json")"
if [ -z "$user_bin_path" ]; then
	echo "set variables"
	exit 1
fi

# Directorys for the production server. Dependencies in JobSubmission.py
if [ ! -d "$production_dir" ]; then
	mkdir "$production_dir"
fi
if [ ! -d "$production_dir/logs" ]; then
	mkdir "$production_dir/logs"
fi
if [ ! -d "$production_dir/results" ]; then
	mkdir "$production_dir/results"
fi
if [ ! -d "$production_dir/job_submissions" ]; then
	mkdir "$production_dir/job_submissions"
fi

# Python
if [ ! -d "$python_env_path" ]; then
	mkdir "$python_env_path"
fi

python3 -m venv "$python_env_path"

# python3 -m pip install --upgrade pip

# shellcheck disable=SC1091,SC1090
source "$python_env_path/bin/activate"

pip install Flask
pip install validate_email
pip install biopython
pip install coolname
pip install ete3
pip install APScheduler
pip install Flask-Limiter
pip install requests-futures
pip install filelock
pip install redis

# R
./system_setup/build_R.sh "$R_bin"

"$R_bin/R/bin/Rscript" ./system_setup/install_packages.R

# mview to build html alignment
wget -nv -P "$user_bin_path/" https://sourceforge.net/projects/bio-mview/files/bio-mview/mview-1.67/mview-1.67.tar.gz

tar -xzf "$user_bin_path/mview-1.67.tar.gz" -C "$user_bin_path"

cd "$user_bin_path/mview-1.67"

./install.pl "$user_bin_path"
