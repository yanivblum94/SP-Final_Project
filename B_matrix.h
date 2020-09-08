/*
 * B_matrix.h
 *
 *  Created on: 21 баев„ 2020
 *      Author: irist
 */

#ifndef B_MATRIX_H_
#define B_MATRIX_H_

/*spmat* create_B(spmat* A, int* ranks, int m, int size);*/

void calc_Bg(spmat* A, spmat* Bg, int* g, int size, int* ranks, int m);

double calc_fi(spmat* Bg, int i);

void calc_Bg_bar(spmat* Bg, int size);

double calc_B_eigen_pair(spmat* B, double *eigenvector, int n);

#endif /* B_MATRIX_H_ */
