/*
 * module_alg.h
 *
 *  Created on: 18 ×‘×�×•×’×³ 2020
 *      Author: irist
 */
#include "spmat.h"

#ifndef MODULE_ALG_H_
#define MODULE_ALG_H_

typedef struct _list_of_lists{

	/*node*/
	 node node;

	/*next pointer*/
	struct _list_of_lists* next;

}list_of_lists;

typedef struct _node{

	/*value*/
	 int val;

	/*next pointer*/
	struct _node* next;

}node;

double calc_Q(int* s, spmat* B, int n);

int mult_vectors(const double *v1, const double *v2, int n);

int* divition_to_2(spmat* B_g_bar, int n);

double calculate_deltaQ(int* s, spmat* B);

void indices_start(int* indices,int n);

int unmoved_start(int* unmoved,int n,int* s);

void modularity_maximization(spmat* BgHat , int* s);




#endif /* MODULE_ALG_H_ */
