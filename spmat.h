

#ifndef SPMAT_H_
#define SPMAT_H_

#ifndef _SPMAT_H
#define _SPMAT_H

typedef struct _spmat {
	/* Matrix size (n*n) */
	int		n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void	(*add_row)(struct _spmat *A, const double *row, int i);

	/* Frees all resources used by A */
	void	(*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void	(*mult)(const struct _spmat *A, const double *v, double *result);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void	*private;
} spmat;

typedef struct _linked_list{

	/*value*/
	double val;

	/*column index*/
	int col;

	/*next pointer*/
	struct _linked_list* next;

}linked_list;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate(int n);

spmat* create_B(spmat* A, int* ranks, int m, int size);

#endif

#endif /* SPMAT_H_ */
