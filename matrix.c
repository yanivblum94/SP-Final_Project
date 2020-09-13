/*
 * matrix.c
 *
 *      Author: irist
 */
#include "spmat.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>

matrix* allocate_matrix(spmat* A, int size, int* k, double* km, int* g){
	matrix* matrixB = (matrix*)malloc(sizeof(matrix));
	if(matrixB == NULL){
		free(matrixB);
		return NULL;
	}
	matrixB->A = A;
	matrixB->size = size;
	matrixB->c = 0.0;
	matrixB->k = k;
	matrixB->km = km;
	matrixB->g = g;

	return matrixB;
}

void free_matrix(matrix* Matrix){
	Matrix->A->free(Matrix->A);
	free(Matrix->g);
	free(Matrix->k);
	free(Matrix->km);
	free(Matrix);
}

double* calc_f(const struct _matrix* B){
	int j,i;
	double sum;
	double* f;
	f = (double*)calloc(B->size, sizeof(double));
	for(j = 0; j < B->size; j++){
		sum = 0.0;
		if(B->g[j] != 0){
			for(i = 0; i< B->size; i++){
				if(B->g[i] != 0){
				sum += B->km[j] * B->k[i];
			}
			if(((linked_list**)(B->A->private))[j] != NULL){
				linked_list *curr;
				curr = ((linked_list**)(B->A->private))[j];
				while(curr != NULL){
					if(B->g[curr->col] != 0){
					sum += curr->val;
					}
					curr = curr->next;
				}
				free(curr);
			}
			}
	}
		f[j] = sum;
	}
	return f;
}

void mult_vector_with_Kmatrix(const struct _matrix* B, const int* v, double* result){
	int i, j;
	double dotproduct;
	for(i = 0; i < B->size; i++){
		dotproduct = 0.0;
		if(B->g[i]!= 0){
			for(j=0; j < B->size; j++){
				if(B->g[j] != 0){
					dotproduct += B->km[j] * B->k[i] * v[j];
		}
			}
				}
	}
	result[i] = dotproduct;
}

void mult_vector_with_I(const struct _matrix* B, double* result){
	int i;
	for(i = 0; i < B->size; i++){
		if(B->g[i]!=0){
		result[i] = B->c;
		}
	}
}

void sum_3_vectors(const double* v1, const double* v2,double* f,  double* result, int size){
	int i;
	for(i = 0; i < size; i++){
		result[i] = v1[i] - v2[i] - f[i];
	}
}

void mult_vector_with_sparse(const struct _matrix* B, const int* v, double* result){
	int i, j;
	double dotproduct;
	linked_list *currlist;
	for(j = 0 ; j < B->size; j++){
		dotproduct = 0.0;
		if(B->g[j]!=0){
			for(i = 0; i < B->size; i++){
				currlist = ((linked_list**)(B->A -> private))[i];
				while(currlist != NULL && currlist->col <= j){
					if(currlist->col == i  &&  B->g[currlist->col]!=0 ){
						dotproduct += (double)(currlist->val * v[i]);
					}
					currlist = currlist->next;
				}
				free(currlist);
			}
		}
		result[j] = dotproduct;
	}
}

void mult_vector_with_matrix(const struct _matrix* B, const int* s, double* result){
	double *v1, *v2, *f;
	f = calc_f(B);
	v1 = (double*)calloc( B->size,sizeof(double));
	v2 = (double*)calloc( B->size,sizeof(double));
	mult_vector_with_sparse(B, s, v1);
	mult_vector_with_Kmatrix(B, s, v2);
	sum_3_vectors(v1, v2, f, result, B->size);
	free(v1);
	free(v2);
	free(f);
}

double mult_vectors_double(const double* v1, const double* v2, int n){
	int i;
	double result = 0.0;
	for(i = 0; i < n; i++){
		result += (*(v1+i))*(*(v2+i));
	}
	return result;
}

double mult_vectors_int(const double* v1, const int* v2, int n){
	int i;
	double result = 0.0;
	for(i = 0; i < n; i++){
		result += (*(v1+i))*(*(v2+i));
	}
	return result;
}

void mult_sparse_with_vector(const struct _matrix* B, const double* v, double* result){
	int i,j;
	double sum;
	linked_list *curr;
	for(i = 0; i < B->size; i++){
		sum=0.0;
		if(B->g[i] != 0){
			curr = ((linked_list**)(B->A->private))[i];
			while(curr != NULL){
				j = curr->col;
				if(B->g[j] != 0){
					sum += curr->val * v[j];
				}
				curr = curr->next;
			}
				free(curr);
			}
		result[i] = sum;
	}
}

void mult_Kmatrix_with_vector(const struct _matrix* B, const double* v, double* result){
	int i, j;
	double dotproduct;
	for(i = 0; i < B->size; i++){
		dotproduct = 0.0;
		if(B->g[i]!= 0){
			for(j=0; j < B->size; j++){
				if(B->g[j] != 0){
					dotproduct += B->km[i] * B->k[j] * v[j];
		}
			}
				}
	}
	result[i] = dotproduct;
}



void sum_4_vectors(const double* v1, const double* v2,const double* v3, const double* f,  double* result, int size){
	int i;
	for(i = 0; i < size; i++){
		result[i] = v1[i] - v2[i] - f[i] + v3[i];
	}
}

void mult_shifted_matrix_with_vector(const struct _matrix* B, const double* v, double* result){
	double *v1, *v2, *v3, *f;
	f = calc_f(B);
	v1 = (double*)calloc( B->size,sizeof(double));
	v2 = (double*)calloc( B->size,sizeof(double));
	v3 = (double*)calloc( B->size,sizeof(double));
	mult_sparse_with_vector(B, v, v1);
	mult_Kmatrix_with_vector(B, v, v2);
	mult_vector_with_I(B, v3);
	sum_4_vectors(v1, v2, v3, f, result, B->size);
	free(v1);
	free(v2);
	free(v3);
	free(f);
}

double calc_norm_1(const struct _matrix* Matrix){
	int i, j, size, flag;
	linked_list *currlist;
	double max, colsum;
	max = 0.0;
	colsum = 0.0;
	size = Matrix->size;
	for(j = 0; j < size; j++){
		colsum = 0.0;
		for(i = 0; i < size; i++){
			flag = 0;
				currlist =  ((linked_list**)(Matrix->A -> private))[i];
				while((currlist != NULL) && (currlist->col <= j)){
					if(currlist->col == j){
						colsum += fabs(currlist->val - ((Matrix->km[i])*(double)(Matrix->k[j])));
						flag = 1;
					}
					currlist = currlist->next;

				}
				if(!flag){
					colsum += fabs((Matrix->km[i])*(Matrix->k[j]));
				}
		}
		if(colsum > max){
			max = colsum;
		}
	}
	return max;
}
