#!/bin/bash

##
## This shell script needs the following programs:
##   program.c
##   monte_carlo.R
##
## In order to generate raw data again, please type
##   rm -f raw_data.tex
##   make raw_data.tex
##

S1=11
S2=12
S3=13

simplex (){
    echo "Numerical experiments for simplex probabilities." > raw_data1.txt
    for i in $(seq 2 1 $1)
    do 
	echo $i
	echo dim=$i >> raw_data1.txt
	/usr/bin/time -o tmp -v ./a.out $S1 $i  > tmp2 
	cat tmp  >> raw_data1.txt
	cat tmp2 >> raw_data1.txt
	echo "" >> raw_data1.txt
    done
    rm tmp tmp2
}

simplicial_cone(){
    file="raw_data2.txt"
    echo "Numerical experiments for normal probabilities of simplicial cone." \
    > $file
    for i in $(seq 2 1 $1)
    do 
	echo $i
	echo dim=$i >> $file
	/usr/bin/time -o tmp -v ./a.out $S2 $i  > tmp2 
	cat tmp  >> $file
	cat tmp2 >> $file
	echo "" >> $file
    done
    rm tmp tmp2
}

simplex_small (){
    echo "Numerical experiments for simplex probabilities(small case)." > raw_data3.txt
    for i in $(seq 2 1 $1)
    do 
	echo $i
	echo dim=$i >> raw_data3.txt
	/usr/bin/time -o tmp -v ./a.out $S3 $i  > tmp2 
	cat tmp  >> raw_data3.txt
	cat tmp2 >> raw_data3.txt
	echo "" >> raw_data3.txt
    done
    rm tmp tmp2
}


##
## main
##

# maximal dimension is 10.
MAXDIM=10

simplex $MAXDIM
simplicial_cone $MAXDIM
simplex_small $MAXDIM

R -f monte_carlo.R