#!/bin/bash

# ./read.sh p0 n n_grp bic iter

p0=$1
n=$2
ngrp=$3
size=$(($ngrp - 1))
bic=$4
TOTAL=$(($p0*$size))

./gather.sh $1 $2 $3 $4 $5

# echo "Here1"

### Begin Truth ###

tLEN=3

for i in $(seq 0 $(($(($p0 * $size)) - 1)))
do
    trueB[$i]=0
done

if [ $ngrp -eq 2 ]
then
    trueB[0]=-1.2

    trueB[4]=0.5

    trueB[8]=1

    trueGRP=(0 4 8)
elif [ $ngrp -eq 3 ]
then
    trueB[0]=-1.2
    trueB[1]=1.8

    trueB[4]=0.5
    trueB[5]=1

    trueB[8]=1
    trueB[9]=1

    trueGRP=(0 2 4)
elif [ $ngrp -eq 4 ]
then
    trueB[0]=-1.2
    trueB[1]=1.8
    trueB[2]=0.25
 
    trueB[12]=0.5
    trueB[13]=-0.5
    trueB[14]=1

    trueB[69]=0.5
    trueB[70]=1
    trueB[71]=-1

    trueGRP=(0 4 23)
elif [ $ngrp -eq 5 ]
then
    trueB[0]=-1.2
    trueB[1]=1.8
    trueB[2]=0.25
    trueB[3]=-1.1

    trueB[12]=0.5
    trueB[13]=-0.5
    trueB[14]=1
    trueB[15]=0.25

    trueB[72]=0.5
    trueB[73]=1
    trueB[74]=-1
    trueB[75]=0.25

    trueGRP=(0 3 18)
fi

### End Truth ###
# echo "Here2"

for a in $(seq 0 1)
do
    filename="results.$1.$2.$3.$a.$4.out"
    cnt=0
    numGRP=0
    N=0
    ss=0
    avgT=0
    avgCov=0
    avgCor0=0
    avgInc0=0
    avgExt=0
    avgMSE=0
    avgSize=0
# echo "Here3"

    while read line
    do
        IFS=',' read -a arrIN <<< "$line"
        if [ $(($cnt % 3)) -eq 0 ]
        then
# echo "Here4 $cnt"
            eLEN=${arrIN[0]}
            Size[$N]=${arrIN[0]}
            for i in $(seq 0 $(($eLEN - 1)))
            do
                estGRP[$i]=${arrIN[$i+1]}
            done

            for i in $(seq 0 $(($tLEN - 1)))
            do 
                for j in $(seq 0 $(($eLEN - 1)))
                do 
                    if [ ${trueGRP[$i]} = ${estGRP[$j]} ]
                    then
                        numGRP=$(($numGRP + 1))
                    fi
                done
            done
            if [ $numGRP -eq $tLEN ]
            then
                Cov[$N]=1
                Inc0[$N]=0
                if [ $eLEN -eq $tLEN ]
                then
                    Ext[$N]=1
                    Cor0[$N]=1
                else
                    Ext[$N]=0
                    Cor0[$N]=$(echo "scale=9; (($p0 - $tLEN) - ($eLEN - $tLEN)) / ($p0 - $tLEN)" | bc)
                fi
            else
                Inc0[$N]=$(echo "scale=9; ($tLEN - $numGRP) / $tLEN" | bc)
                Cov[$N]=0
                Ext[$N]=0 
                Cor0[$N]=$(echo "scale=9; (($p0 - $tLEN) - ($p0 - ($eLEN - $numGRP))) / ($p0 - $tLEN)" | bc)
            fi
            numGRP=0
            avgSize=$(echo "scale=9; $avgSize + ${Size[$N]}" | bc)
            avgCov=$(echo "scale=9; $avgCov + ${Cov[$N]}" | bc)
            avgCor0=$(echo "scale=9; $avgCor0 + ${Cor0[$N]}" | bc)
            avgInc0=$(echo "scale=9; $avgInc0 + ${Inc0[$N]}" | bc)
            avgExt=$(echo "scale=9; $avgExt + ${Ext[$N]}" | bc)
        elif [ $(($cnt % 3)) -eq 1 ]
        then
# echo "Here5 $cnt"
            MSE=0
            for i in $(seq 0 $(($(($p0 * $size)) - 1)))
            do
# echo "${arrIN[$i]} ${trueB[$i]} $i"
                num=$(echo "scale=9; ${arrIN[$i]} - ${trueB[$i]}" | bc)
                num=$(echo "scale=9; $num^2" | bc)
                MSE=$(echo "scale=9; $MSE + $num" | bc)
                ((ss++))
            done
            arrMSE[$N]=$(echo "scale=9; $MSE/$ss" | bc)
            ss=0
            MSE=0
            avgMSE=$(echo "scale=9; $avgMSE + ${arrMSE[$N]}" | bc)
        else
# echo "Here6 $cnt"
            T[$N]=${arrIN[0]}
            avgT=$(echo "scale=9; $avgT + ${T[$N]}" | bc)
            ((N++))
        fi   
        ((cnt++))
    done < $filename

    avgSize=$(echo "scale=9; $avgSize/$N" | bc)
    avgCov=$(echo "scale=9; $avgCov/$N" | bc)
    avgCor0=$(echo "scale=9; $avgCor0/$N" | bc)
    avgInc0=$(echo "scale=9; $avgInc0/$N" | bc)
    avgExt=$(echo "scale=9; $avgExt/$N" | bc)
    avgT=$(echo "scale=9; $avgT/$N" | bc)
    avgMSE=$(echo "scale=9; $avgMSE/$N" | bc)

    sdT=0
    sdCov=0
    sdCor0=0
    sdInc0=0
    sdExt=0
    sdMSE=0
    sdSize=0

    tmp1=0
    tmp2=0
    tmp3=0
    tmp4=0
    tmp5=0
    tmp6=0
    tmp7=0

    for i in $(seq 0 $(($N - 1)))
    do
        tmp1=$(echo "scale=9; ${Cov[$i]} - $avgCov" | bc)
        sdCov=$(echo "scale=9; $sdCov + $tmp1^2" | bc)
        tmp2=$(echo "scale=9; ${Cor0[$i]} - $avgCor0" | bc)
        sdCor0=$(echo "scale=9; $sdCor0 + $tmp2^2" | bc)
        tmp3=$(echo "scale=9; ${Inc0[$i]} - $avgInc0" | bc)
        sdInc0=$(echo "scale=9; $sdInc0 + $tmp3^1" | bc)
        tmp4=$(echo "scale=9; ${Ext[$i]} - $avgExt" | bc)
        sdExt=$(echo "scale=9; $sdExt + $tmp4^2" | bc)
        tmp5=$(echo "scale=9; ${T[$i]} - $avgT" | bc)
        sdT=$(echo "scale=9; $sdT + $tmp5^2" | bc)
        tmp6=$(echo "scale=9; ${arrMSE[$i]} - $avgMSE" | bc)
        sdMSE=$(echo "scale=9; $sdMSE + $tmp6^2" | bc)
        tmp7=$(echo "scale=9; ${Size[$i]} - $avgSize" | bc)
        sdSize=$(echo "scale=9; $sdSize + $tmp7^2" | bc)
        tmp1=0
        tmp2=0
        tmp3=0
        tmp4=0
        tmp5=0
        tmp6=0
        tmp7=0
    done

    sdCov=$(echo "scale=9; $sdCov/$N" | bc)
    sdCor0=$(echo "scale=9; $sdCor0/$N" | bc)
    sdInc0=$(echo "scale=9; $sdInc0/$N" | bc)
    sdExt=$(echo "scale=9; $sdExt/$N" | bc)
    sdT=$(echo "scale=9; $sdT/$N" | bc)
    sdMSE=$(echo "scale=9; $sdMSE/$N" | bc)
    sdSize=$(echo "scale=9; $sdSize/$N" | bc)

    sdCov=$(echo "scale=9; sqrt($sdCov)" | bc)
    sdCor0=$(echo "scale=9; sqrt($sdCor0)" | bc)
    sdInc0=$(echo "scale=9; sqrt($sdInc0)" | bc)
    sdExt=$(echo "scale=9; sqrt($sdExt)" | bc)
    sdT=$(echo "scale=9; sqrt($sdT)" | bc)
    sdMSE=$(echo "scale=9; sqrt($sdMSE)" | bc)
    sdSize=$(echo "scale=9; sqrt($sdSize)" | bc)

    if [ $bic -eq 1 ]
    then
        echo "BIC = Small BIC calculation" > "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    else
        echo "BIC = Large BIC calculation" > "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    fi

    echo "Total number of Groups = $p0" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "Total Number of Parameters = $TOTAL" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "Group Size=$ngrp" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"

    if [ $a -eq 0 ]
    then
        echo "Variance Structure = AR1" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    else
        echo "Variance Structure = IID" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    fi

    echo "Coverage probability = $avgCov ($sdCov)"  >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "Percentage of correct zeros = $avgCor0 ($sdCor0)" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "Percentage of incorrect zeros = $avgInc0 ($sdInc0)" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "Exact select probability = $avgExt ($sdExt)" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "Model Size = $avgSize ($sdSize)" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "MSE = $avgMSE ($sdMSE)" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"
    echo "Total time = $avgT ($sdT)" >> "summary_result.$p0.$n.$ngrp.$a.$bic.out"

done

exit 0
