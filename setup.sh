#!/usr/bin/env bash

if module 2>/dev/null; then
	OLDMODULEPATH=$MODULEPATH

	MODULEPATH="$MODULEPATH:/home/skulkarni/privatemodules"
	module load gcc/5.1.0 # seems to prevent a very strange bug
	module load samtools
	module load htslib/git-head
	module load yaml
	module load gengetopt

	MODULEPATH=$OLDMODULEPATH

	echo "Modules loaded."
fi
