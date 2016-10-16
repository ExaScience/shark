#!/bin/bash

export MPICH_ASYNC_PROGRESS=1

if [ $# -ne 4 ]        
then                   
	echo Usage: ./run.sh nb_processes nb_threads size max_iterations 1>&2
	exit 1
fi

tol=1e-15
maxit=$4
t=$2
n=$3
ppn=$1
r=1

rm timing*.dat

for SIZE in $n #250 500 1000
do
    for Pn in {1..1..1}
    do

	head -n $Pn $PBS_NODEFILE > nodefile

	 MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=noverbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r -c pipe_chronopoulos_cg 2>&1 
	
	 MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=noverbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r -c pipe_cg 2>&1
	
	 MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=noverbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r -c pipe_gmres -l 1 2>&1 

	 MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=noverbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r -c pipe_gmres -l 2 2>&1 

	 MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=noverbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r -c pipe_gmres -l 3 2>&1 

	 MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=noverbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r -c pipe_gmres -l 4 2>&1	

	 MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=noverbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn $ppn $HOME/exa2ct/shark/build/main -n $SIZE -m $maxit -e $tol -t $t -r $r -c pipe_gmres -l 5 2>&1 
	
	done

done


