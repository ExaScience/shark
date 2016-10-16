#!/bin/bash

MV2_ENABLE_AFFINITY=0 GOMP_CPU_AFFINITY=0,1,2,3,6,7,8,9 KMP_AFFINITY=verbose,granularity=fine,proclist=[0,1,2,3,6,7,8,9] mpirun -f $PBS_NODEFILE -ppn 1 ./cg
