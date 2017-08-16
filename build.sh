#!/usr/bin/env bash

# die on any error
set -e

source ./setup.sh

init=
debug=
makeflags=

badopt() {
	echo "Bad option: $1" >&2
}

while true
do
	case $1 in
		--init | -i ) init=1 ; shift  ;;
		--debug | -d ) debug=1 ; shift ;;
		--clean ) clean=1 ; shift ;;
		"" ) break ;;
		-- ) shift ; break ;;
		* ) badopt $1; exit 1 ;;
	esac
done

if [[ $init ]]
then
	autoreconf --install
	cd build
	../configure
else
	cd build
fi

if [[ $debug ]]
then
	makeflags="$makeflags CFLAGS=\"-O0 -DDEBUG\""
fi

if [[ $clean ]]
then
	makeflags="$makeflags clean"
fi

eval "make $makeflags"
