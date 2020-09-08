/*
 * B_matrix.c
 *
 *  Created on: 18 баев„ 2020
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

/*spmat* create_B(spmat* A, int* ranks, int m, int size){
	spmat* B;
	double* row;
	int i,j,col;
	B = spmat_allocate(size);
	for(i=0; i<size; i++){
		linked_list *curr = ((linked_list**)(A->private))[i];
		row = (double*)calloc(size, sizeof(double));
		col = curr->col;
		for(j=0; j< size; j++){
			if(j==col){
				row[j] = 1 - ((ranks[j]*ranks[i])/m);
				curr = curr->next;
			}
			else{
				row[j] = 0 - ((ranks[j]*ranks[i])/m);
			}
		}
		B->add_row(B, row, i);
		free(row);
	}
	return B;
}*/


/*void copy_matrix(spmat* A, spmat* B, int n){
	int i;
	double *row;
	linked_list *curr;
	B = spmat_allocate(n);
	for(i = 0; i < n; i++){
		curr = ((linked_list**)(A->private))[i];
		row = (double*)malloc(n*sizeof(double));
		while(curr != NULL){
			row[curr->col] = curr->val;
		}
		B->add_row(B, row, i);
		free(row);
	}
}*/

void calc_Bg(spmat* A, spmat* Bg, int* g, int size, int* ranks, int m){
	int i, col;
	linked_list *curr;
	double *row;
	for(i = 0; i < size; i++){
		if(g[i] != 0){
			curr = ((linked_list**)(A->private))[i];
			row = (double*)calloc(size, sizeof(double));
			while(curr != NULL){
				col = curr->col;
				if(g[col] != 0){
					row[col] = curr->val - ((ranks[col]*ranks[i])/m);
				}
				curr = curr->next;
			}
			Bg->add_row(Bg, row, i);
			free(row);
		}
	}
	free(curr);
}


double calc_fi(spmat* Bg, int i){
	double sum;
	linked_list *list = ((linked_list**)(Bg->private))[i];
	sum = 0;
	while(list != NULL){
		sum += list->val;
		list = list->next;
	}
	return sum;
}

void calc_Bg_bar(spmat* Bg, int size){
	int i;
	linked_list *curr;
	for(i = 0; i < size; i++){
		int fi = calc_fi(Bg, i);
		curr = ((linked_list**)(Bg->private))[i];
		while((curr != NULL) && (curr->col <= i)){
			if(curr->col == i){
				curr->val = curr->val - fi;
			}
			curr = curr->next;
		}
	}
}


double calc_B_eigen_pair(spmat* B, double *eigenvector, int n){
	double *currvector;
	double eigen_val;
	spmat *shifted_matrix = (spmat*)malloc(sizeof(spmat));
	currvector = (double*)calloc(n, sizeof(double));
	initialize_random_vector(currvector, n);
	calc_shift(B, shifted_matrix);
	calc_eigen(shifted_matrix, currvector, eigenvector, n);
	eigen_val = calc_eigen_val(shifted_matrix, eigenvector, n);
	eigen_val -= calc_norm_1(B);
	free_in_list(shifted_matrix);
	free(currvector);
	return eigen_val;

}



