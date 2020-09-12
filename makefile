FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: matrix.o cluster.o eigen_pair.o module_alg.o spmat.o
	gcc -g matrix.o cluster.o eigen_pair.o module_alg.o spmat.o -o cluster $(LIBS)

matrix.o: matrix.c matrix.h
	gcc $(FLAGS) -c matrix.c
cluster.o: cluster.c
	gcc $(FLAGS)  -c cluster.c 

eigen_pair.o: eigen_pair.c
	gcc $(FLAGS) -c eigen_pair.c 

module_alg.o: module_alg.c
	gcc $(FLAGS) -c module_alg.c 

spmat.o: spmat.c
	gcc $(FLAGS) -c spmat.c 


clean:
	rm -f -rf *.o prog





