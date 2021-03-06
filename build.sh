#!/usr/bin/env bash

# die on any error
set -e

source ./setup.sh

init=
debug=
scan=
makeflags=
makepreamble=
cflags=

badopt() {
	echo "Bad option: $1" >&2
}

while true
do
	case $1 in
		--init | -i ) init=1 ; shift  ;;
		--debug | -d ) debug=1 ; shift ;;
		--clean ) clean=1 ; shift ;;
		--scan ) scan=1 ; shift ;;
		"" ) break ;;
		-- ) shift ; break ;;
		* ) badopt $1; exit 1 ;;
	esac
done

if [[ $init ]]
then
	autoreconf --install
	./configure
fi

if [[ $debug ]]
then
	cflags="-g -Og -DDEBUG -DHASH_DEBUG=1"
else
	cflags="-O3 -DNDEBUG"
fi

makeflags="$makeflags CFLAGS=\"$cflags\""

if [[ $clean ]]
then
	makeflags="$makeflags clean"
fi

if [[ $scan ]]
then
	makepreamble="scan-build "
fi

eval "${makepreamble}make $makeflags"
