CC = gcc
CFLAGS = -Wall -g -I./
LDFLAGS = -L./
ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

default: libwrap_tetinter.so

all: libwrap_tetinter.so

libwrap_tetinter.so:
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -Wl,-rpath,$(ROOT_DIR) -o libwrap_tetinter.so -fPIC -fopenmp plane_normal_tetrahedron_intersect.c

call_tetinter: libwrap_tetinter.so
	$(CC) $(CFLAGS) $(LDFLAGS)  call_libwrap_tetinter.c -lwrap_tetinter -lm -Wl,-rpath,${PWD} -o call_tetinter

test: call_tetinter
	./call_tetinter

clean:
	rm -f libwrap_tetinter.so call_tetinter
