/*
 * matrix.h
 *
 *      Author: irist
 */

#ifndef MATRIX_H_
#define MATRIX_H_

typedef struct _matrix{
	spmat* A;

	int size;

	double c;

	int* k;

	double* km;

	int* g;
}matrix;

matrix* allocate_matrix(spmat* A, int size, int* k, double* km, int* g);

void free_matrix(matrix* Matrix);

double* calc_f(const struct _matrix* B);

void mult_vector_with_Kmatrix(const struct _matrix* B, const int* v, double* result);

void mult_vector_with_I(const struct _matrix* B, double* result);

void sum_3_vectors(const double* v1, const double* v2,double* f,  double* result, int size);

void mult_vector_with_sparse(const struct _matrix* B, const int* v, double* result);

void mult_vector_with_matrix(const struct _matrix* B, const int* s, double* result);

double mult_vectors_double(const double* v1, const double* v2, int n);

double mult_vectors_int(const double* v1, const int* v2, int n);

void mult_sparse_with_vector(const struct _matrix* B, const double* v, double* result);
void mult_Kmatrix_with_vector(const struct _matrix* B, const double* v, double* result);

void sum_4_vectors(double* v1, double* v2,double* v3, double* f, double* result, int size);

void mult_shifted_matrix_with_vector(const struct _matrix* B, const double* v, double* result);

double calc_norm_1(const struct _matrix* Matrix);


#endif /* MATRIX_H_ */
