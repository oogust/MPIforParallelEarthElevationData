CC = mpicc
LDLIBS = -lm -pg
CFLAGS = -g -O3
CFLAGS += -Wall

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi) 
LDLIBS += $(shell pkg-config --libs scalapack-openmpi) 

ALL: model validate

model: model.o harmonics.o
validate: validate.o harmonics.o 
model.o: harmonics.h model.h
quality.o: harmonics.h
harmonics.o: harmonics.h

install_scalapack : 
	sudo-g5k apt-get install libscalapack-openmpi-dev

CR : CR.tex
	pdflatex CR.tex

.PHONY: clean

clean:
	rm -f CR.aux CR.log CR.out model validate *.o