# Makefile for SPECFEMX
# REVISION
#   HNG, Jan 17,2014

default: createdir specfem-x

all: createdir default partmesh specfem-x pspecfem-x

createdir:
	(mkdir -p bin; mkdir -p input; mkdir -p output; mkdir -p partition; mkdir -p tmp)

clean:
	(cd src; make clean)

debug:
	(cd src; make $@)

test:
	(cd src; make $@)

specfem-x: 
	(cd src; make $@)

pspecfem-x: 
	(cd src; make $@)

partmesh: 
	(cd src; make $@)

