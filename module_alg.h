/*
 * module_alg.h
 *
 *      Author: irist
 */
#include "spmat.h"

#ifndef MODULE_ALG_H_
#define MODULE_ALG_H_

typedef struct _node{

	/*value*/
	 int val;

	/*next pointer*/
	struct _node* next;

}node;

typedef struct _list_of_lists{

	/*node*/
	 node* node;

	/*next pointer*/
	struct _list_of_lists* next;

}list_of_lists;



double calc_Q(int* s, spmat* B, int n);

double mult_vectors_double(const double *v1, const double *v2, int n);

double mult_vectors_int(const double *v1, const int *v2, int n);

int* divition_to_2(spmat* B_g_bar, int n);

double calculate_deltaQ(int* s, spmat* B);

void indices_start(int* indices,int n);

int unmoved_start(int* unmoved,int n,int* s);

void modularity_maximization(spmat* BgHat , int* s);

list_of_lists* divide_network(spmat* B, int n);




#endif /* MODULE_ALG_H_ */
