#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>


spmat* read_mat(FILE *input, int* ranks, int size){
	int i,k,temp,j;
	double* row;
	spmat* A;
	A = spmat_allocate(size);
	for(i=0; i<size; i++){
		fread(&k, sizeof(int), 1, input);
		ranks[i] = k;
		row = calloc(k, sizeof(double));
		fread(&row, sizeof(int), k, input);
		A->add_row(row);
	}
	return A;
}

int calc_M(int* ranks,int size){
	int i, m=0;
	for(i=0; i< size; i++){
		m+= ranks[i];
	}
	return m;
}


int main(int argc, char* argv[]){
	FILE *input;
	spmat* A,B;
	int* ranks;
	int size,m;
	input = fopen(argv[1]);
	fread(&size, sizeof(int), 1, input);
	ranks = calloc(size, sizeof(int));
	A = read_mat(input,ranks,size);
	m = calc_M(ranks,size);
	B = create_B(A, ranks, m, size);
}
