
#!/bin/bash

while read line  
do  

str=`echo $line | cut -f1 -d ":" `
if [ "$str" == 'n' ]; then
size=`echo $line | cut -f2 -d ":" `
elif [ "$str" == 'maxit' ]; then
iteration=`echo $line | cut -f2 -d ":" `
elif [ "$str" == 'nprocs' ]; then 
nprocs=`echo $line | cut -f2 -d ":" `
elif [ "$str" == 'nthrds' ]; then
nthrds=`echo $line | cut -f2 -d ":" `
elif [ "$str" == 'err' ]; then
error=`echo $line | cut -f2 -d ":" `
elif [ "$str" == 'runtime' ]; then
runtime=`echo $line | cut -f2 -d ":" `
echo $nprocs,$nthrds,$size,$iteration,$error,$runtime >>parsed_output.dat
fi;

done < $1
