/*
 * cluster.c
 *
 *  Created on: 17 באוג׳ 2020
 *      Author: irist
 */
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
	int n;
	double* row;
	spmat* A;
	A = spmat_allocate(size);
	for(i = 0; i < size; i++){
		n = fread(&k, sizeof(int), 1, input);
		if(n != 1){
			printf("Couldn't open the file");
		}
		ranks[i] = k;
		row = (double*)calloc(k, sizeof(double));
		n = fread(&row, sizeof(int), k, input);
		if(n != k){
			printf("Error in reading a row");
		}
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
	spmat *A, *B;
	int *ranks;
	int size,m, n;
	input = fopen(argv[1], "r");
	if(input == NULL){
		printf("Error in reading the file");
		return 1;
	}
	n = fread(&size, sizeof(int), 1, input);
	if(n != 1){
		printf("Error in reading the size of the matrix");
	}
	ranks = (int*)calloc(size, sizeof(int));
	if(ranks == NULL){
		printf("Error in allocation of memory");
		return 1;
	}
	A = read_mat(input,ranks,size);
	fclose(input);
	m = calc_M(ranks,size);
	B = create_B(A, ranks, m, size);
	
	free(ranks);
	free_in_list(A);
	free_in_list(B);
	return 0;
}

