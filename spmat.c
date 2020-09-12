/*
 * spmat.c
 *
 *  Created on: 17 баев„ 2020
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

void forceStop(const char* func,const int line){
	printf("\n ~ forced to stop in %s, line %d", func, line);
	exit(EXIT_SUCCESS);
}
/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate(int n){
	spmat* matrix = (spmat*)malloc(sizeof(spmat));
	if(matrix==NULL){
			matrix->free(matrix);
			return NULL;
	}
	matrix->private = (linked_list**)malloc(n*sizeof(linked_list*));
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
	return matrix;

}


void add_row_in_list(struct _spmat *A, const double *row, int i){
	int j,flag;
	linked_list *elem, *curr;
	flag =1;
	((linked_list**)(A->private))[i] = NULL;
	for(j = 0; j < A->n; j++){
		if(row[j] != 0.0){
			elem = (linked_list*)malloc(sizeof(linked_list));
			if(elem==NULL){
				A->free(A);
				return;
			}
			elem->val = row[j];
			elem->col = j;
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
}


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





