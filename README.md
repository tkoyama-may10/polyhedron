## Synopsis

This program computes the probability content of simplex regions
and simplicial cone regions with a multivariate normal distribution
by the holonomic gradient method (HGM).
For details, see `note-en.tex`.

## Code Example

In /bin/bash/, the following command will compute the probability content
with a multivariate normal distribution
of the simplex given by x1 > 0, x2 > 0, and -x1-x2+10>0.

`./a.out 1 2  1 0 0  0 1 0 -1 -1 10`

Here, the first argument switches the behavior of the program.
The second argument is the dimension of the simplex
and the remaining arguments are the coefficients of the inequalities.


## Motivation

On this project, we release programs and raw data used in our paper [1].

## Installation

Download the zip file `simplex-master.zip`.
Then run the following commands:
```
unzip simplex-master.zip
cd simplex-master/
gcc program.c -lm -lgsl -lblas 
```
For compiling `program.c`, you need the GNU scientific library, BLAS and LAPACK.


## API Reference

1. T.Koyama, 
``
Holonomic gradient method for the probability content of a simplex region
with a multivariate normal distribution
'', 
http://arxiv.org/abs/1512.06564, 2015.
2. GNU scientific library, http://www.gnu.org/software/gsl/

## Tests

```
$ ./a.out 1 2  1 0 0  0 1 0  -1 -1 10
d=2
...
Probability = 0.25
$ ./a.out 2 2 1 0 0 -1 1 0 
d=2
...
Probability = 0.125
```


## Contributors

Tamio Koyama

## License

GPL

