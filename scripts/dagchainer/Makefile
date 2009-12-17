
PROG=dagchainer
PYTHON_INCLUDE=/usr/include/python2.6

#CCFLAGS=-Wno-deprecated
#CCFLAGS=-O3 -s -fomit-frame-pointer -msse2 -Wno-deprecated
CCFLAGS=-O3 -s -fomit-frame-pointer -Wno-deprecated
all: 
	g++ -o ${PROG} dagchainer.cpp ${CCFLAGS}

python:
	cython cdagline.pyx && \
		gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I${PYTHON_INCLUDE} -c cdagline.c -o cdagline.o && \
		gcc -pthread -shared -Wl,-O2 -Wl,-Bsymbolic-functions cdagline.o -o cdagline.so && \
		rm cdagline.o && rm cdagline.c 

