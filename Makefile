all: a.out note-en.pdf

a.out: program.c
	gcc program.c -lm -lgsl -lblas -llapack #-O0 -Wall

note-en.pdf:note-en.tex
	pdflatex note-en.tex && pdflatex note-en.tex

