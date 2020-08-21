/*
 * module_alg.h
 *
 *  Created on: 18 באוג׳ 2020
 *      Author: irist
 */
#include "spmat.h"

#ifndef MODULE_ALG_H_
#define MODULE_ALG_H_

double calc_Q(int* s, spmat* B, int n);

int mult_vectors(const double *v1, const double *v2, int n);

int* divition_to_2(spmat* B_g_bar, int n);




#endif /* MODULE_ALG_H_ */
