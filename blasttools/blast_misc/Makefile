
all:
	cython blast_misc.pyx
	gcc -O3 -I/usr/include/python2.5 -c -o blast_misc.o blast_misc.c
	gcc -shared blast_misc.o -o blast_misc.so


#clean:
#    rm -f blast_misc.c 
#	rm -rf blast_misc.o
#	rm -rf blast_misc.so
