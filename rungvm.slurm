#!/usr/bin/env bash

#SBATCH --job-name=gvm
#SBATCH --array=1-24
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --output=output/%A_%a.out

CHRS=(DUMMY 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
echo "Chromosome $CHR"

cd $SLURM_SUBMIT_DIR

source ./setup.sh
module load yaml
module load samtools
module load htslib

GVM_EXEC=$PWD/build/gvm

SC="/scratch/skulkarni"
PERL5LIB="$PERL5LIB:/home/rhalperin/perl5lib/lib"
CONF=/home/skulkarni/programs/gvm/conf.yaml

JOBNUM=$SLURM_ARRAY_JOB_ID

SCDIR=$SC/$JOBNUM

echo "GVM: $GVM_EXEC"

mkdir $SCDIR
rm $SC/latest
ln -s $SCDIR $SC/latest

sed "s/\\\$OUTFILE/$JOBNUM\/out/g" $CONF > $SCDIR/conf.yaml

touch $SCDIR/still_running_chr$CHR
perf stat $GVM_EXEC $SCDIR/conf.yaml $CHR
rm $SCDIR/still_running_chr$CHR


