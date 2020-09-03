/*
 * cluster.c
 *
 *      Author: irist
 */
#include "spmat.h"
#include "module_alg.h"
#include "eigen_pair.h"
#include "B_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>


spmat* read_mat(FILE *input, int* ranks, int size){
	int i,k,j;
	int n;
	double* row;
	int * row_tmp;
	spmat* A;
	A = spmat_allocate(size);
	printf("entered read_mat \n");
	for(i = 0; i < size; i++){
		n = fread(&k, sizeof(int), 1, input);
		if(n != 1){
			printf("Couldn't open the file");
		}
		ranks[i] = k;
		row = (double*)calloc(k, sizeof(double));
		row_tmp = (int*)calloc(k, sizeof(int));
		n = fread(&row_tmp, sizeof(int), k, input);
		if(n != k){
			printf("Error in reading a row");
		}
		for(j=0;j<k;j++){
			row[j] = (double)row_tmp[j];
		}
		A->add_row(A, row, i);
		free(row);
		free(row_tmp);
	}
	return A;
}


int calc_M(int* ranks,int size){
	int i, m=0;
	for(i=0; i< size; i++){
		m += ranks[i];
	}
	return m;
}

int calc_set_size(node* set){
	node *curr = set;
	int res =0;
	if(set == NULL){
		return res;
	}
	while(curr != NULL){
		res++;
		curr = curr->next;
	}
	free(curr);
	return res;
}

int calc_num_sets(list_of_lists* sets){
	int res;
	list_of_lists *curr =  sets;
	res=0;
	while(curr != NULL){
		curr = curr->next;
	}
	free(curr);
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

void list_to_arr(int* arr, node* set, int size){
	int i;
	node *curr = set;
	for(i=0; i<size; i++){
		arr[i] = curr->val;
		curr = curr->next;
	}
	qsort(arr,size, sizeof(int), compare);
}


int main(int argc, char* argv[]){
	FILE *input, *output;
	spmat *A, *B;
	int *ranks;
	int size, m, n, sets_num, temp, asserter;
	list_of_lists *sets;
	printf("start running the program \n");
	if(argc != 3){
		printf("More/Less then 2 parameters have been given");
		return 1;
	}

	input = fopen(argv[1], "r");
	if(input == NULL){
		printf("Error in reading the file");
		return 1;
	}
	printf("created input \n");
	n = fread(&size, sizeof(int), 1, input);
	if(n != 1){
		printf("Error in reading the size of the matrix");
	}
	printf("size=  %d \n" , size);
	ranks = (int*)calloc(size, sizeof(int));
	if(ranks == NULL){
		printf("Error in allocation of memory");
		return 1;
	}
	A = read_mat(input,ranks,size);
	printf("read matrix A \n");
	fclose(input);
	m = calc_M(ranks,size);
	B = create_B(A, ranks, m, size);
	sets = divide_network(B,size);
	sets_num = calc_num_sets(sets);
	output = fopen(argv[2], "w");
	temp = fwrite(&sets_num, sizeof(int),1,output);/*write num of sets*/
	if(temp!=1){
		return 1;
	}
	while(sets!=NULL){
		int *set;
		temp = calc_set_size(sets->node);
		asserter = fwrite(&temp, sizeof(int),1,output);
					if(asserter!=1){
						return 1;
					}
		set = (int*)calloc(temp, sizeof(int));
		list_to_arr(set, sets->node, temp);
		asserter = fwrite(&set, sizeof(int),temp,output);
		if(asserter!=temp){
			return 1;
		}
		sets = sets->next;
		free(set);
	}
	free(ranks);
	free(sets);
	fclose(output);
	free_in_list(A);
	free_in_list(B);
	return 0;
}

