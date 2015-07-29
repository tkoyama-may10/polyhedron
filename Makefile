all: a.out 

a.out: program.c
	gcc program.c -lm -lgsl -lblas -llapack #-O0 -Wall


