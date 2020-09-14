
#include "spmat.h"
#include "module_alg.h"
#include "eigen_pair.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>
#include "matrix.h"

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

void poweriteration(matrix *Matrix, double *currvector, double *nextvector, int n){
	int i;
	double norm;
	mult_shifted_matrix_with_vector(Matrix, currvector, nextvector);
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

void calc_eigen(matrix* Matrix, double* currvector, double* nextvector, int n){
	int i;
	poweriteration(Matrix, currvector, nextvector, n);
	while(check(currvector, nextvector, n) != 1){
		for(i = 0; i < n; i++){
			currvector[i] = nextvector[i];
		}
		poweriteration(Matrix, currvector, nextvector, n);
	}
}

double calc_eigen_val(matrix *Matrix, double *eigenvector, int n){
	double eigenval;
	double numerator, denominator;
	double *matrix_by_vec = (double*)calloc(n, sizeof(double));
	mult_shifted_matrix_with_vector(Matrix, eigenvector, matrix_by_vec);
	numerator = mult_vectors_double(matrix_by_vec, eigenvector, n);
	denominator = mult_vectors_double(eigenvector, eigenvector, n);
	eigenval = numerator/denominator;
	return eigenval;

}

