#!/bin/bash

set -e

python_env_backend_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['python_env_backend_path'])" \
	< "./parameters_backend_local.json")"
if [ -z "$python_env_backend_path" ]; then
	echo "set variables"
	exit 1
fi

python_env_seqSel_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['python_env_seqSel_path'])" \
	< "./parameters_backend_local.json")"
if [ -z "$python_env_seqSel_path" ]; then
	echo "set variables"
	exit 1
fi

blast_bin_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['blast_bin_path'])" \
	< "./parameters_backend_local.json")"
if [ -z "$blast_bin_path" ]; then
	echo "set variables"
	exit 1
fi

clustal_bin_path="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['clustal_bin_path'])" \
	< "./parameters_backend_local.json")"
if [ -z "$clustal_bin_path" ]; then
	echo "set variables"
	exit 1
fi

work_dir="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['work_dir'])" \
	< "./parameters_backend_local.json")"
if [ -z "$work_dir" ]; then
	echo "set variables"
	exit 1
fi

if [ ! -d "$work_dir" ]; then
	mkdir "$work_dir"
fi

# python env seqSel
if [ ! -d "$python_env_seqSel_path" ]; then
	mkdir "$python_env_seqSel_path"
fi

python3 -m venv "$python_env_seqSel_path"

# shellcheck disable=SC1091,SC1090
source "$python_env_seqSel_path/bin/activate"

pip install biopython
pip install ete3
pip install six

pip freeze > ./templates/python_requirements.txt

deactivate

# python env backend
if [ ! -d "$python_env_backend_path" ]; then
	mkdir "$python_env_backend_path"
fi

python3 -m venv "$python_env_backend_path"

# shellcheck disable=SC1091,SC1090
source "$python_env_backend_path/bin/activate"

pip install jinja2
pip install biopython
pip install ete3
pip install six

deactivate

# blast
if [ ! -d "$blast_bin_path" ]; then
	mkdir "$blast_bin_path"
fi

blast_version=$(curl https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 2> /dev/null | \
    grep -o 'ncbi-blast-.*x64-linux\.tar\.gz"' | \
    awk -F "-" '{print $3}')

wget -nv -P "$blast_bin_path" "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${blast_version}-x64-linux.tar.gz"

tar -xzf "$blast_bin_path/ncbi-blast-$blast_version-x64-linux.tar.gz" -C "$blast_bin_path"

cp "$blast_bin_path/ncbi-blast-$blast_version/bin/"* "$blast_bin_path"

rm -rf "$blast_bin_path/ncbi-blast-$blast_version"
rm "$blast_bin_path/ncbi-blast-$blast_version-x64-linux.tar.gz"

PATH="$PATH:$blast_bin_path"

# shellcheck disable=SC1091
source ./system_setup/update_blast_db.sh

# clustal

wget -nv -P "$clustal_bin_path" \
	http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64
