#!/bin/bash

##
## This shell script needs the following programs:
##   program.c
##   monte_carlo.R
##   dist_func.R
##
## In order to generate raw data again, please type
##   gcc program.c -lm -lgsl -lblas -llapack
##   ./numerical_experiment.sh
##

S1=11
S2=12
S3=13
S4=14

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

dist_func (){
echo \# Raw data for graphs > tmp

for dim  in 2 6 10
do
    echo \#dim= $dim >> tmp
    for c in 0.0 0.5 1.0 1.5  2.0 2.5  3.0 3.5  4.0  4.5 5.0 5.5  6.0 
    do
	./a.out $S4 $dim $c >> tmp
    done
done
cat tmp > raw_data4.txt
rm tmp
R -f dist_func.R
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

dist_func