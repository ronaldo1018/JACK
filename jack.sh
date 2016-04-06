#!/bin/sh
clear
begin_time=`date '+%s'`
# Show OpenCL Platforms & Devices
./jack q

# Number of Traces or Particles
TRACESNUM=(1000 10000 1000000)
# Collect grid resolution
COLLECT=(100 200)
# Trace resolution
MAXSTEP=(0.01 0.02 0.05 0.1)
# Incident particle energy
ENERGY=(200 150 100 50)

# for CPU code
NTHREAD=(1 2 4 8)
THREADSIZE=${#NTHREAD[*]}

# for OpenCL code
ITEM_PER_WG=(0 0 0)
ITEM_PER_BATCH=(1000 10000 250000)

echo
printf "How many OpenCL Devices will be run : "
read -n 1 OCLNUM

if [[ $OCLNUM = "" ]];
then
    OCLNUM=0
fi
echo
for (( i = 0; i < $OCLNUM; i++ ))
do
    printf "Select OpenCL platform : "
    read word 
    PLAT[$i]=$word
    printf "Select OpenCL Device : "
    read word 
    DEV[$i]=$word
done

clear

echo "Combination of traces are  ${TRACESNUM[@]}"
echo "Combination of collect are ${COLLECT[@]}"
echo "Combination of steps are   ${MAXSTEP[@]}"
echo "Combination of energy are  ${ENERGY[@]}"
echo
echo "Combination of CPU threads are  ${NTHREAD[@]}"
echo
printf "Select OpenCL device is(are) "

if [ $OCLNUM -eq 0 ]
then    printf " none."
fi

echo

for (( i = 0; i < $OCLNUM; i++ ))
do
    echo "    Platform : ${PLAT[$i]}  Device : ${DEV[$i]}"
done

echo
echo "Press ENTER to get start...."
read tmp

if [ -f ./jack_script.csv ]
then     rm ./jack_script.csv
         echo "Delete jack_script.csv ...."
fi  
    
    
###### OpenCL Devices #####    
for (( a = 0; a < $OCLNUM; a++ ))
do
    for (( i = 0; i < ${#TRACESNUM[@]}; i++ ))
    do
        for (( j = 0; j < ${#COLLECT[@]}; j++ ))
        do
            for (( m = 0; m < ${#MAXSTEP[@]}; m++ ))
            do
                for (( n = 0; n < ${#ENERGY[@]}; n++ ))
                do
                    echo
                    echo "========================================================"
                    echo "Running OpenCL Device => Platform:${PLAT[$a]}  Device:${DEV[$a]}"
                    echo "   Nmc=${TRACESNUM[$i]} ; dC=${COLLECT[$j]} ; dS=${MAXSTEP[$m]} ; E=${ENERGY[$n]}"
                    echo "   ./jack so ${TRACESNUM[$i]} ${PLAT[$a]} ${DEV[$a]} ${ITEM_PER_WG[$i]} ${ITEM_PER_BATCH[$i]} ${COLLECT[$j]} ${MAXSTEP[$m]} ${ENERGY[$n]}"
                    ./jack so ${TRACESNUM[$i]} ${PLAT[$a]} ${DEV[$a]} ${ITEM_PER_WG[$i]} ${ITEM_PER_BATCH[$i]} ${COLLECT[$j]} ${MAXSTEP[$m]} ${ENERGY[$n]}
                done
            done   
        done
    done
done 

###### CPU core #####    
for (( a = 0; a < $THREADSIZE; a++ ))
do
    for (( i = 0; i < ${#TRACESNUM[@]}; i++ ))
    do
        for (( j = 0; j < ${#COLLECT[@]}; j++ ))
        do
            for (( m = 0; m < ${#MAXSTEP[@]}; m++ ))
            do
                for (( n = 0; n < ${#ENERGY[@]}; n++ ))
                do
                    echo
                    echo "========================================================"
                    echo "Running CPU code => Thread(s) = ${NTHREAD[$a]} "
                    echo "   Nmc=${TRACESNUM[$i]} ; dC=${COLLECT[$j]} ; dS=${MAXSTEP[$m]} ; E=${ENERGY[$n]}"
                    echo "   ./jack sc ${TRACESNUM[$i]} ${NTHREAD[$a]} ${COLLECT[$j]} ${MAXSTEP[$m]} ${ENERGY[$n]}"
                    ./jack sc ${TRACESNUM[$i]} ${NTHREAD[$a]} ${COLLECT[$j]} ${MAXSTEP[$m]} ${ENERGY[$n]}
                done
            done   
        done
    done
done

end_time=`date '+%s'`
total_time=$end_time-$begin_time
ss=$((total_time%60))
mm=$(((total_time/60)%60))
hh=$((total_time/3600))
echo "========================================================"
echo "Script Finished"
printf "Total Time = %i:%02i:%02i" $hh $mm $ss
echo
echo "Results saved : jack_script.csv"
echo "========================================================"
echo
   
