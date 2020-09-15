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

/**
 * reads the input file into a sparse matrix A
 * @param A - the matrix in which we add the rows
 * @param input - the input graph file
 * @param ranks - the vector f ranks which we change in the method
 * @param size - size if matrix
 */
void read_mat(spmat* A, FILE *input, int* ranks, int size){
	int i,k;
	int n;
	int *row_tmp;
	/*char *func = "red_mat";*/
	printf("entered read_mat \n");
	for(i = 0; i < size; i++){
		n = fread(&k, sizeof(int), 1, input);
		if(n != 1){
			printf("Couldn't read the %d", i);
		}
		ranks[i] = k;
		/*printf("ranks[%d] = %d\n" , i,ranks[i] );*/
		if(k > 0){
			row_tmp = (int*)calloc(k, sizeof(int));
			if(row_tmp == NULL){
				printf("Couldn't allocate memory");
			}
			n = fread(row_tmp, sizeof(int), k, input);
			if(n != k){
				printf("Error in reading a row");
			}
		A->add_row(A, row_tmp, i, k);
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
			printf("  val:  %d" ,  currlist->val);
			currlist = currlist->next;
		}
		printf("\n");
	}

}

/**
 * calculates the value M
 * @param ranks - the ranks vector of the nodes
 * @param size - size of ranks
 * @return
 */
int calc_M(int* ranks,int size){
	int i, m=0;
	for(i=0; i< size; i++){
		m += ranks[i];
	}
	printf("m = %d", m);
	return m;
}
/**
 * calculates size of a divided set in order to write in the output
 * @param set - the referance to the first node of the "set list"
 * @return the number of nodes in the set
 */
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
/**
 * calculates how many sets we have in the division
 * @param sets - the list of lists which we have to count the num of lists
 * @return - number of sets
 */
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

/**
 * helps the sorting
 * @param a - int a
 * @param b - int b
 * @return 0 for equal, 1 if a>b, -1 else
 */
int compare( const void* a, const void* b)
{
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}
/**
 * takes a list and makes it into an array
 * @param arr - the array we want
 * @param set - the list we get
 * @param size - size of array
 */
void list_to_arr(int* arr, node* set, int size){
	int i;
	node *curr = set;
	for(i=0; i<size; i++){
		arr[i] = curr->val;
		curr = curr->next;
	}
	qsort(arr,size, sizeof(int), compare);
}
/**
 * Calculates the value of Ki/M
 * @param ranks - vector of ranks
 * @param ranks_m - the result vector that will hold Ki/M in each ith cell
 * @param m - M
 * @param size - size of ranks
 */
void calc_ranks_m(int* ranks, double* ranks_m, int m, int size){
	int i;
	for(i = 0; i < size; i++){
		ranks_m[i] = ((double)(ranks[i]))/((double)m);
		printf("rank[%d] = %f" , i,ranks_m[i] );
	}
}


int main(int argc, char* argv[]){
	FILE *input, *output;
	spmat *A;
	int *ranks;
	double *ranks_m;
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
	A = spmat_allocate(size);
	read_mat(A,input,ranks,size);
	print_mat(A);
	printf("read matrix A \n");
	printf("matrix A size:  %d \n", A->n);
	fclose(input);
	/*print_mat(A);*/
	ranks_m = (double*)calloc(size, sizeof(double));
	m = calc_M(ranks,size);
	calc_ranks_m(ranks, ranks_m, m, size);
	sets = divide_network(A, size, ranks, ranks_m);
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
		asserter = fwrite(set, sizeof(int),temp,output);
		if(asserter!=temp){
			return 1;
		}
		sets = sets->next;
		free(set);
	}
	free(ranks);
	free(ranks_m);
	free(sets);
	fclose(output);
	free_in_list(A);
	return 0;
}

