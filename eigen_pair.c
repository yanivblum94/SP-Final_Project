
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

/**
 * initializes random vector
 * @param currvector - the vector we give values to
 * @param n - size of vector
 */
void initialize_random_vector(double *currvector, int n){
	int i;
	srand(time(NULL));
	for(i = 0; i < n; ++i){
		currvector[i] = (double)rand();
	}
}

/**
 * calculates the norm of a vector
 * @param nextvector - the vector we calculate the norm to
 * @param n - size of vector
 * @return the norm of the vector
 */
double calcnorm(double *nextvector, int n){
	int i;
	double sum;
	sum = 0.0;
	for(i = 0; i < n; i++){
		sum += pow(nextvector[i],2.0);
	}
	return sqrt(sum);
}

/**
 * PI that gets the vector closer to the eigen vector of the Matrix
 * @param Matrix - the matrix we want the eigenvector of
 * @param currvector - the current eigen vector we have
 * @param nextvector - the vector after the iteration
 * @param n - size od matrix
 */
void poweriteration(matrix *Matrix, double *currvector, double *nextvector, int n){
	int i;
	double norm;
	mult_shifted_matrix_with_vector(Matrix, currvector, nextvector);
	norm = calcnorm(nextvector, n);
	for(i = 0; i < n; i++){
		nextvector[i] = (double) nextvector[i]/norm;
	}
}

/**
 * Boolean function checking if the current eigen vec we have is "close enpugh" according to the epsilon defined to the next one we received from PI
 * @param currvector - the current vector
 * @param nextvector - the next vector we got
 * @param n - size of vector
 * @return - 1 if we close enough, 0 else
 */
int check(double *currvector, double *nextvector, int n){
	int i;
	for(i = 0; i < n; ++i){
		if(fabs((double)(currvector[i] - nextvector[i])) >= EPSILON){
			return 0;
		}
	}
	return 1;
}


/**
 * calculates the eigen vector
 * @param Matrix - the matrix we calculate the eigen to
 * @param currvector - the current vector we have
 * @param nextvector - the next vector we get from PI - the final vector
 * @param n - size of vector
 */
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
/**
 * calculates the leading eigen value
 * @param Matrix - the matrix we calculate the eigen val to
 * @param eigenvector - the eigen vector
 * @param n - size f eigen vector
 * @return - eigen value
 */
double calc_eigen_val(matrix *Matrix, double *eigenvector, int n){
	double eigenval;
	double numerator, denominator;
	double *matrix_by_vec = (double*)calloc(n, sizeof(double));
	mult_shifted_matrix_with_vector(Matrix, eigenvector, matrix_by_vec);
	numerator = mult_vectors_double(matrix_by_vec, eigenvector, n);
	denominator = mult_vectors_double(eigenvector, eigenvector, n);
	eigenval = numerator/denominator;
	free(matrix_by_vec);
	return eigenval;

}

