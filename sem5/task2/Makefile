CC=gcc
CFLAGS=-O3 -march=native
SWIG=swig
PYTHONPATH=/usr
PYTHONVERSION=3.6

all: _mymodule.so

qr.o: qr.c qr.h
	$(CC) qr.c -c -fPIC $(CFLAGS)

mymodule_wrap.c: mymodule.i qr.h
	$(SWIG) -python mymodule.i

mymodule.py: mymodule.i qr.h
	$(SWIG) -python mymodule.i

mymodule_wrap.o: mymodule_wrap.c
	$(CC) mymodule_wrap.c -c -I$(PYTHONPATH)/include/python$(PYTHONVERSION)/ -I/home/azamat42/.local/lib/python$(PYTHONVERSION)/site-packages/numpy/core/include -fPIC $(CFLAGS)

_mymodule.so: qr.o mymodule_wrap.o
	$(CC) mymodule_wrap.o qr.o -ldl -lm -o _mymodule.so -shared $(CFLAGS)

.PHONY: clean

clean:
	rm -f qr.o mymodule_wrap.o _mymodule.so mymodule.py mymodule_wrap.c
