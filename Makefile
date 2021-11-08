IDIR =./iohb1.0

default: all

all: bicg
	mpicc gradiente_biconjugado.c iohb1.0/iohb.o -Wno-pointer-to-int-cast -Wno-format-security -I$(IDIR) -o gradiente_biconjugado
bicg:
	mpicc -Wno-pointer-to-int-cast -Wno-format-security -I$(IDIR) -c $(IDIR)/iohb.c
