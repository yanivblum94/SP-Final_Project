/*
 * cluster.c
 *
 *  Created on: 17 ×‘×�×•×’×³ 2020
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

int calc_set_size(linked_list* set){
	linked_list curr = set;
	int res =0;
	if(set == NULL){
		return res;
	}
	while(set!= NULL){
		res++;
		curr = curr.next;
	}
	return res;
}

int calc_num_sets(linked_list** sets, int size){
	int i,res;
	res=0;
	for(i=0; i<size; i++){
		if(sets[i] !=NULL){
			res++;
		}
	}
	return res;
}

int compare( const void* a, const void* b)
{
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

void list_to_arr(int* arr, linked_list* set, int size){
	int i;
	linked_list curr = set;
	for(i=0; i<size; i++){
		arr[i] = curr.val;
		curr = curr.next;
	}
	qsort(arr,size, sizeof(int), compare);
}


int main(int argc, char* argv[]){
	FILE *input, *output;
	spmat *A, *B;
	int *ranks;
	int size,m, n,sets_num,temp,i,asserter;
	linked_list** sets = (linked_list**)malloc(n*sizeof(linked_list*));
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

	
	
	sets_num = calc_num_sets(sets,size);
	output = fopen(argv[2], "w");
	temp = fwrite(&sets_num, sizeof(int),1,output);/*write num of sets*/
	if(temp!=1){
		return 1;
	}
	for(i=0; i<size; i++){
		temp = calc_set_size(sets[i]);
		if(temp>0){
			int* set = (int*)calloc(temp, sizeof(int));
			list_to_arr(set, sets[i], temp);
			asserter = fwrite(&set, sizeof(int),temp,output);
			if(asserter!=temp){
				return 1;
			}
		}
	}
	free(ranks);
	fclose(output);
	free_in_list(A);
	free_in_list(B);
	return 0;
}

