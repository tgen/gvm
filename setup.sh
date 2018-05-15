#!/usr/bin/env bash

if module 2>/dev/null; then
	OLDMODULEPATH=$MODULEPATH

	MODULEPATH="$MODULEPATH:/home/skulkarni/privatemodules"
	module load gcc/6.1.0 # seems to prevent a very strange bug
	module load samtools
	module load htslib/git-head
	module load yaml
	module load gengetopt
	module load gsl/2.4-sid
        module load python/3.6.0

	# TEMPORARY STOPGAP MEASURE (HPC ticket 39196)
	export LIBRARY_PATH="/packages/gsl/2.4/lib:$LIBRARY_PATH"

	MODULEPATH=$OLDMODULEPATH

	echo "Modules loaded."
fi
