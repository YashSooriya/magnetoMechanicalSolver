#! /bin/bash

NS=( 180 90 45 23 )
F1=( 10 20 40 80 )
F2=( 50 100 200 400 )
for index in "${!NS[@]}"
do
    matlab -nodisplay -nodesktop -r "run_parallel(${F1[index]},${F2[index]}, ""PODI"");exit;">& output_PODI_${NS[index]}_parallel.txt
    echo "Done N_s = ${NS[index]}"
done
