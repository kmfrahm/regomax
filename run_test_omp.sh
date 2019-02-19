echo "$1 0.15 10 10000 100 $2 $3" > INPUTtestreduce_omp\_$1
time (./net_testgreduced_omp.out < INPUTtestreduce_omp\_$1 > OUTPUTtestreduce_omp\_$1 ) 2> timetestreduce_omp\_$1
