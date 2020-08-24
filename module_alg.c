/*
 * module_alg.c
 *
 *  Created on: 18 באוג׳ 2020
 *      Author: irist
 */

#include "spmat.h"
#include "eigen_pair.h"
#include "B_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>
#include "module_alg.h"


int mult_vectors(const double *v1, const double *v2, int n){
	int i;
	double result = 0.0;
	for(i = 0; i < n; i++){
		result += (*(v1+i))*(*(v2+i));
	}
	return result;
}

double calc_Q(int* s, spmat* B, int n){
	double q;
	double *result = (double*)malloc(n*sizeof(double));
	mult_vector_with_matrix(B, s, result);
	q = 0.5*mult_vectors(result,s);
	free(result);
	return q;
}

int division_to_2(spmat* B, int* g, int n){
	double eigen_val;
	double *eigenvector;
	int *s, i, is_divisible;
	is_divisible = 1;
	eigenvector = (double*)malloc(n*sizeof(double));
	s = (int*)malloc(n*sizeof(int));
	spmat *Bg = (spmat*)malloc(sizeof(spmat));
	calc_Bg(B, Bg, g, n);
	calc_Bg_bar(Bg, n);
	eigen_val = calc_B_eigen_pair(Bg, eigenvector, n);
	if(eigen_val <= 0){//The group is indivisible
		is_divisible = 0;
	}
	else{
		for(i = 0; i < n; i++){
			if(eigenvector[i] > 0){
				s[i] = 1;
			}
			else{
				s[i] = -1;
			}
		}
	}
	if(calc_Q(s, B, n) <= 0){//The group is indivisible
		is_divisible = 0;
	}
	else{
		for(i = 0; i < n; i++){
			if(g[i] != 0){
				g[i] = s[i];
			}
		}
	}
	
	
	free(eigenvector);
	free(s);
	free_in_list(Bg);
	return is_divisible;
}


