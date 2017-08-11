#!/usr/bin/env bash

# die on any error
set -e

source ./setup.sh

if [[ "$1" != "--rebuild" ]]
then
	autoreconf --install
	cd build
	../configure
else
	cd build
fi

make
