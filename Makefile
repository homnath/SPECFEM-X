# Makefile for SPECFEMX
# REVISION
#   HNG, Jan 17,2014

default: createdir specfemx

all: createdir default partmesh specfemx pspecfemx

createdir:
	(mkdir -p bin; mkdir -p input; mkdir -p output; mkdir -p partition; mkdir -p tmp)

clean:
	(cd src; make clean)

debug:
	(cd src; make $@)

test:
	(cd src; make $@)

specfemx: 
	(cd src; make $@)

pspecfemx: 
	(cd src; make $@)

partmesh: 
	(cd src; make $@)

