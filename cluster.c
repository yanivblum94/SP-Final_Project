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

void forceStop(const char* func,const int line){
	printf("\n ~ forced to stop in %s, line %d", func, line);
	exit(EXIT_SUCCESS);
}

spmat* read_mat(FILE *input, int* ranks, int size){
	int i,k,j;
	int n;
	double *row;
	int *row_tmp;
	spmat *A;
	char *func = "red_mat";
	A = spmat_allocate(size);
	printf("entered read_mat \n");
	for(i = 0; i < size; i++){
		n = fread(&k, sizeof(int), 1, input);
		if(n != 1){
			printf("Couldn't read the %d", i);
		}
		printf("k of ith row %d %d /n" , k, i);
		ranks[i] = k;
		if(k>0){
		row = (double*)calloc(k, sizeof(double));
		if(row == NULL){
			printf("Couldn't allocate memory");
			forceStop(func, 41);
		}
		row_tmp = (int*)calloc(k, sizeof(int));
		if(row_tmp == NULL){
			printf("Couldn't allocate memory");
			forceStop(func, 46);
		}
		n = fread(row_tmp, sizeof(int), k, input);
		if(n != k){
			printf("Error in reading a row");
		}
		if(i == 1){
			/*forceStop(func, 53);*/
		}
		for(j = 0; j < k; j++){
			/*forceStop(func, 57);*/
			printf("row[j] = %f", row[j]);
			/*forceStop(func, 59);*/
			printf("row_tmp[j] = %d", row_tmp[j]);
			/*forceStop(func, 59);*/
			/*if(j == 0){
				forceStop(func, 54);
			}*/
			/*forceStop(func, 59);*/
			row[j] = (double) row_tmp[j];
			/*forceStop(func, 57);*/
		}
		/*forceStop(func, 47);*/
		A->add_row(A, row, i);
		free(row);
		free(row_tmp);
	}
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
	spmat *A;
	int *ranks;
	int size, m, n, sets_num, temp, asserter;
	list_of_lists *sets;
	char *func = "main";
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
	/*forceStop(func, 138);*/
	printf("read matrix A \n");
	fclose(input);
	m = calc_M(ranks,size);
	/*forceStop(func, 142);
	forceStop(func, 144);*/
	sets = divide_network(A, size, m, ranks);
	forceStop(func, 146);
	sets_num = calc_num_sets(sets);
	forceStop(func, 148);
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
	return 0;
}

