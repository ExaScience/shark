#! /bin/bash

function sweep()
{
    lib=$1
    runcmd=$2
    heatcmd=$3
    for n in 1 2 4
    do
        for t in 1 # 2 4 8 16 32
        do
            set -x
            $runcmd -n $n $heatcmd -t $t >$lib.$n.$t.log 2>&1
            set +x 
        done
    done
}

N=8192
I=20
sweep gpi "$HOME/exa2ct/gpi2-release/bin/gaspi_run -m $HOME/nodes" "./heat.gpi -n $N -i $I"
sweep mpi "mpirun -f $HOME/nodes" "./heat.mpi -n $N -i $I"
