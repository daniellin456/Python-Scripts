CC=gcc
CFLAGS=-std=c11
FC=gfortran
FFLAGS=


all:
	$(CC) $(CFLAGS) -lm Create_Table.c Export_EoS_Table.c -o ../bin/Export_EoS_Table
#	$(FC) $(FFLAGS) readtable.f95 -o ../bin/readtable
clean:
	rm -f ../bin/Export_EoS_Table
#	rm -r ../bin/readtable
