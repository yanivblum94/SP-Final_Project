
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>


/**
 *  Allocates a new linked-lists sparse matrix of size n
 * @param n - size of matrix
 * @return spamt with allocated memory
 */
spmat* spmat_allocate(int n){
	spmat* matrix = (spmat*)calloc(1, sizeof(spmat));
	if(matrix==NULL){
			matrix->free(matrix);
			return NULL;
	}
	matrix->private = (linked_list**)calloc(n, sizeof(linked_list*));
	if(matrix->private==NULL){
				matrix->free(matrix);
				return NULL;
		}
	matrix->add_row = &add_row_in_list;
	if(matrix->add_row==NULL){
				matrix->free(matrix);
				return NULL;
		}
	matrix->free = &free_in_list;
	if(matrix->free==NULL){
				matrix->free(matrix);
				return NULL;
		}
	matrix->mult = &mult_matrix_with_double_vector;
	if(matrix->mult==NULL){
				matrix->free(matrix);
				return NULL;
		}
	matrix->n = n;
	return matrix;

}

/**
 * adds row to the matrix based on the input file described for the project
 * @param A - the spmat he rows are added to
 * @param row - array of numbers to be added
 * @param i - the row number
 * @param k - size of array
 */
void add_row_in_list(struct _spmat *A, const int *row, int i, int k){
	int j,flag;
	linked_list *elem, *curr;
	flag =1;
	((linked_list**)(A->private))[i] = NULL;
	for(j = 0; j < k; j++){
		elem = (linked_list*)calloc(1, sizeof(linked_list));
		if(elem==NULL){
			A->free(A);
			return;
		}
		elem->val = 1;
		elem->col = row[j];
		if(flag){
			curr = elem;
			((linked_list**)(A->private))[i] = curr;
			flag=0;
		}
		else{
				curr->next = elem;
				curr= curr->next;
				curr->next = NULL;
		}
	}
}

/**
 * frees the spmat
 * @param A - the mat that is freed
 */
void free_in_list(struct _spmat *A){
	int i;
	linked_list *currlist, *temp;
	linked_list **l = (linked_list**)(A -> private);
	if(l==NULL){
			A->free(A);
			return;
	}
	for(i = 0; i < A->n; i++){
		currlist = l[i];
		while(currlist != NULL){
			temp = currlist;
			currlist = currlist->next;
			free(temp);
		}
	}
	free(A->private);
	free(A);
}

/**
 * performs mult of matrix with vector of double values
 * @param A - The matrix being multiplied
 * @param v - the vector being multiplied
 * @param result - the result vector which is changed in the method
 */
void mult_matrix_with_double_vector(const struct _spmat *A, const double *v, double *result){
	int i, j, n;
	double dotproduct;
	linked_list *currlist;
	n = A->n;
	for(i = 0 ; i < n; i++){
		dotproduct = 0.0;
		currlist = ((linked_list**)(A -> private))[i];
		while(currlist != NULL){
			j = currlist->col;
			dotproduct += (double)(currlist->val * v[j]);
			currlist = currlist->next;
		}
		result[i] = dotproduct;
	}
}

/**
 * performs mult of matrix with vector of int values
 * @param A - The matrix being multiplied
 * @param v - the vector being multiplied
 * @param result - the result vector which is changed in the method
 */
void mult_matrix_with_int_vector(const struct _spmat *A, const int *v, double *result){
	int i, j, n;
	double dotproduct;
	linked_list *currlist;
	n = A->n;
	for(i = 0 ; i < n; i++){
		dotproduct = 0.0;
		currlist = ((linked_list**)(A -> private))[i];
		while(currlist != NULL){
			j = currlist->col;
			dotproduct += (double)(currlist->val * v[j]);
			currlist = currlist->next;
		}
		result[i] = dotproduct;
	}
}

/**
 * calculates the norm of a matrix
 * @param A - the matrix whose norm is calculated
 * @return the norm
 */

double calc_norm_1_A(const struct _spmat *A){
	int i,j,n;
	linked_list *currlist;
	double max, colsum;
	max = 0.0;
	colsum = 0;
	n = A->n;
	for(j = 0; j < n; j++){
		colsum = 0.0;
		for(i = 0; i < n; i++){
			currlist =  ((linked_list**)(A -> private))[i];
			while((currlist != NULL) && (currlist->col <= j)){
				if(currlist->col == j){
					if((currlist->val) < 0){
						colsum += -(currlist->val);
					}
					else{
						colsum += (currlist->val);
					}
				}
				currlist = currlist->next;
			}
		}
		if(colsum > max){
			max = colsum;
		}
	}
	return max;
}






