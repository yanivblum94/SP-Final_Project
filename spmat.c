/*
 * spmat.c
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
	matrix->mult = &mult_matrix_with_vector;
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


void mult_matrix_with_vector(const struct _spmat *A, const double *v, double *result){
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

void mult_vector_with_matrix(const struct _spmat *A, const double *v, double *result){

}

/*claculate 1-norm of the matrix*/
double calc_norm_1(const struct _spmat *A){
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
					colsum += abs(currlist->val);
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

/*create I*||C|| matrix*/
void create_IC(spmat* i_matrix, int size, double c){
	int i;
	double *row;
	row = (double*)malloc(size*sizeof(double));
	for(i = 0; i < size; i++){
		row[i] = c;
	}
	for(i = 0; i < size; i++){
		add_row(i_matrix, row, size);
	}
}

/*Calculate the shifted matrix C'*/
void calc_shift(const struct _spmat *A, spmat *shifted_matrix){
	int i, n;
	double c = calc_norm_1(A);
	double *row;
	linked_list *currlist1, *currlist2;
	n = A->n;
	shifted_matrix = create_IC(A, n, c);
	row = (double*)malloc(n*sizeof(double));
	for(i = 0; i < n; i++){
		currlist1 = (double*)(shifted_matrix->private)[i];
		currlist2 = (double*)(A->private)[i];
		while(currlist1 != NULL){
			currlist1->val = currlist1->val + currlist2->val;
		}
	}
}







