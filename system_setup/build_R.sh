#!/bin/bash

set -e

R_path="$1"

if [ ! -d "$R_path" ]; then
	mkdir "$R_path"
fi

cd "$R_path"

wget "http://cran.rstudio.com/src/base/R-4/R-4.1.1.tar.gz"
tar xvf "R-4.1.1.tar.gz"
cd R-4.1.1
./configure --prefix="$R_path/R" --with-readline=no --with-x=no
make && make install
