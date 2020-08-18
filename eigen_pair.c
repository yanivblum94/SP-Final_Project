#include "spmat.h"
#include "module_alg.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>

static const double EPSILON = 0.00001;

void initialize_random_vector(double *currvector, int n){
	int i;
	srand(time(NULL));
	for(i = 0; i < n; ++i){
		currvector[i] = (double)rand();
	}
}

double calcnorm(double *nextvector, int n){
	int i;
	double sum;
	sum = 0.0;
	for(i = 0; i < n; i++){
		sum += pow(nextvector[i],2.0);
	}
	return sqrt(sum);
}

void poweriteration(spmat *matrix, double *currvector, double *nextvector, int n){
	int i;
	double norm;
	matrix->mult(matrix, currvector, nextvector);
	norm = calcnorm(nextvector, n);
	for(i = 0; i < n; i++){
		nextvector[i] = (double) nextvector[i]/norm;
	}
}

int check(double *currvector, double *nextvector, int n){
	int i;
	for(i = 0; i < n; ++i){
		if(fabs((double)(currvector[i] - nextvector[i])) >= EPSILON){
			return 0;
		}
	}
	return 1;
}

void calc_eigen(spmat *matrix, double *currvector, double *nextvector, int n){
	int i;
	poweriteration(matrix, currvector, nextvector, n);
	while(check(currvector, nextvector, n) != 1){
		for(i = 0; i < n; i++){
			currvector[i] = nextvector[i];
		}
		poweriteration(matrix, currvector, nextvector, n);
	}
}

double calc_eigen_val(spmat *matrix, double *eigenvector, int n){
	double eigenval;
	double numerator, denominator;
	double *matrix_by_vec;
	matrix_by_vec = mult_matrix_with_vector(matrix, eigenvector, matrix_by_vec);
	numerator = mult_vectors(matrix_by_vec, eigenvector, n);
	denominator = mult_vectors(eigenvector, eigenvector, n);
	eigenval = numerator/denominator;
	return eigenval;

}

