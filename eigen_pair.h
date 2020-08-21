/*
 * eigen_pair.h
 *
 *  Created on: 21 באוג׳ 2020
 *      Author: irist
 */

#ifndef EIGEN_PAIR_H_
#define EIGEN_PAIR_H_

void initialize_random_vector(double *currvector, int n);

double calcnorm(double *nextvector, int n);

void poweriteration(spmat *matrix, double *currvector, double *nextvector, int n);

int check(double *currvector, double *nextvector, int n);

void calc_eigen(spmat *matrix, double *currvector, double *nextvector, int n);

double calc_eigen_val(spmat *matrix, double *eigenvector, int n);



#endif /* EIGEN_PAIR_H_ */
