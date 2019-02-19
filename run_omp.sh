echo "$1 0.15 10 10000 100 $2 $3" > INPUTreduce_omp\_$1\_$2
time (./net_greduced_omp.out < INPUTreduce_omp\_$1\_$2 > OUTPUTreduce_omp\_$1\_$2 ) 2> timereduce_omp\_$1\_$2
