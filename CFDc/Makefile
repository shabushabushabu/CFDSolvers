CC 		= gcc

all: vectorutils.o finitedifference.o fluiddynamics.o main.o
	$(CC) vectorutils.o finitedifference.o fluiddynamics.o main.o -o cavityflow

vectorutils.o: vectorutils.c
	$(CC) -c vectorutils.c -o vectorutils.o

finitedifference.o: finitedifference.c
	$(CC) -c finitedifference.c -o finitedifference.o

fluiddynamics.o: fluiddynamics.c
	$(CC) -c fluiddynamics.c -o fluiddynamics.o

main.o: main.c
	$(CC) -c main.c -o main.o

clean:
	rm *.o cavityflow