/**
 * this module represents the B and Bg matrixes. We divide it into 3 Matrixes - A, KiKj/M and cI
 * KiKj/M is a mult of Ki/M vector and ranks(k) vector
 */
#ifndef MATRIX_H_
#define MATRIX_H_
/**
 * the struct representing the matrix
 */
typedef struct _matrix{
	spmat* A; /*spmat which we read from input*/

	int size; /* size of A (nXn)*/

	double c; /* the norm*/

	int* k; /*ranks vector*/

	double* km; /*Ki/M vector*/

	int* g; /*reprenets Bg - were g[i] is 1 = we count the row and col for Bg*/
}matrix;

matrix* allocate_matrix(spmat* A, int size, int* k, double* km, int* g);

void free_matrix(matrix* Matrix);

void calc_f(const struct _matrix* B, double* f);

void mult_vector_with_I(const struct _matrix *B, const double *v, double *result);

void sum_3_vectors(const double* v1, const double* v2,double* f,  double* result, int size);

double mult_vectors_double(const double* v1, const double* v2, int n);

double mult_vectors_int(const double* v1, const int* v2, int n);

void mult_vector_with_f_int(const struct _matrix* B, const int* v, const double* f, double* result);

void mult_vector_with_f_double(const struct _matrix* B, const double* v, const double* f, double* result);

void mult_sparse_with_vector(const struct _matrix* B, const double* v, double* result);
void mult_Kmatrix_with_vector(const struct _matrix* B, const double* v, double* result);

void sum_4_vectors(const double* v1, const double* v2,const double* v3, const double* f,  double* result, int size);

void mult_shifted_matrix_with_vector(const struct _matrix* B, const double* v, double* result);

double calc_norm_1(const struct _matrix* Matrix);

void mult_matrix_with_vector(const struct _matrix* B, const int* v, double* result);



#endif /* MATRIX_H_ */
