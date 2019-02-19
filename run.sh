echo "$1 0.15 10 10000 100 $2 $3" > INPUTreduce\_$1\_$2
time (./net_greduced.out < INPUTreduce\_$1\_$2 > OUTPUTreduce\_$1\_$2 ) 2> timereduce\_$1\_$2
