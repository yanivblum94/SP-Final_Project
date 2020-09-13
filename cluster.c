/*
 * cluster.c
 *
 *      Author: irist
 */
#include "spmat.h"
#include "module_alg.h"
#include "matrix.h"
/*#include "eigen_pair.h"*/
/*#include "B_matrix.h"*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>


void read_mat(spmat* A, FILE *input, int* ranks, int size){
	int i,k,j;
	int n;
	double *row;
	int *row_tmp;
	char *func = "red_mat";
	printf("entered read_mat \n");
	for(i = 0; i < size; i++){
		n = fread(&k, sizeof(int), 1, input);
		if(n != 1){
			printf("Couldn't read the %d", i);
		}
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
		for(j = 0; j < k; j++){
			row[j] = (double) row_tmp[j];
		}
		A->add_row(A, row, i, k);
		free(row);
		free(row_tmp);
	}
	}
}

void print_mat(spmat* A){
	int i;
	linked_list *currlist;
	printf("entered print_mat \n");
	for(i=0; i<A->n; i++){
		printf("row:  %d ", i);
		currlist = ((linked_list**)(A->private))[i];
		while(currlist != NULL){
			printf("  col:  %d" ,  currlist->col);
			printf("  val:  %f" ,  currlist->val);
			currlist = currlist->next;
		}
		printf("\n");
	}

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

double* calc_ranks_m(int* ranks, int m, int size){
	double *ranks_m;
	int i;
	ranks_m = (double*)malloc(size*sizeof(double));
	for(i = 0; i < size; i++){
		ranks_m[i] = (double)(ranks[i]/m);
	}
	return ranks_m;
}


int main(int argc, char* argv[]){
	FILE *input, *output;
	spmat *A;
	int *ranks;
	double *ranks_m, norm;
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
	A = spmat_allocate(size);
	read_mat(A,input,ranks,size);
	printf("read matrix A \n");
	printf("matrix A size:  %d \n", A->n);
	fclose(input);
	print_mat(A);
	m = calc_M(ranks,size);
	if(m==0){
		printf("m=0 - error!");
		return 1;
	}
	ranks_m = calc_ranks_m(ranks, m, size);
	norm = calc_norm_1_A(A);
	/*forceStop(func, 176);*/
	sets = divide_network(A, size, ranks, ranks_m, norm);
	forceStop(func, 178);
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
		asserter = fwrite(set, sizeof(int),temp,output);
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

