#include "spmat.h"
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
	return q;
}
