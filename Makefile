IDIR =./iohb1.0

default: all

all: bicg
	gcc gradiente_biconjugado.c iohb1.0/iohb.o -Wno-pointer-to-int-cast -Wno-format-security -I$(IDIR) -o gradiente_biconjugado -pg
bicg:
	gcc -Wno-pointer-to-int-cast -Wno-format-security -I$(IDIR) -c $(IDIR)/iohb.c
