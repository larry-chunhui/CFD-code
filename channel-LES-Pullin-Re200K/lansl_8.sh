#!/bin/sh
#___INFO__MARK_BEGIN__
# Welcome to use  EasyCluster V1.6 All Rights Reserved.
#
#___INFO__MARK_END__
#
#$ -S /bin/sh
#$ -N Re200k
#$ -j y
#$ -o ./
#$ -e ./
#$ -cwd
#$ -q fortran.q
#$ -pe mvapi 8-8 

source ~/.bashrc
hash -r
export path=$TMPDIR:$path
/usr/local/mvapi2/bin/mpiexec -launcher rsh -n $NSLOTS -f $TMPDIR/machines ./main

