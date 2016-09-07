#!/bin/bash

export MPICH_ASYNC_PROGRESS=1

if [ $# -ne 3 ]        
then                   
	echo Usage: ./run.sh nb_processes size max_iterations 1>&2
	exit 1
fi

tol=1e-15
maxit=$3
t=
n=$2
ppn=$1
r=1

rm timing*.dat

for SIZE in $n
do
    for Pn in {1..1..1}
    do

	head -n $Pn $PBS_NODEFILE > nodefile

	runcmd="mpirun -hostfile nodefile -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r "
	
	$runcmd -c pipe_chronopoulos_cg 2>&1 
	
	$runcmd -c pipe_cg 2>&1
	
	$runcmd -c pipe_gmres -l 1 2>&1 
	$runcmd -c pipe_gmres -l 2 2>&1 
	$runcmd -c pipe_gmres -l 3 2>&1 
	$runcmd -c pipe_gmres -l 4 2>&1	
	$runcmd -c pipe_gmres -l 5 2>&1 
	
	done

done


