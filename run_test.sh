echo "$1 0.15 10 10000 100 $2 $3" > INPUTtestreduce\_$1
time (./net_testgreduced.out < INPUTtestreduce\_$1 > OUTPUTtestreduce\_$1 ) 2> timetestreduce\_$1
