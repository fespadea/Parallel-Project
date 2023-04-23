all:
	mpicc proj.c -o proj -lm

clean:
	rm *.btr