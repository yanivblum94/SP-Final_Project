/**
 * this module contains all of the algorithms of the network division
 */
#include "spmat.h"
#include "matrix.h"


#ifndef MODULE_ALG_H_
#define MODULE_ALG_H_

/*the struct node represents the head of a list*/
typedef struct _node{

	/*value*/
	 int val;

	/*next pointer*/
	struct _node* next;

}node;

/*this struct represents the division groups that we have*/
typedef struct _list_of_lists{

	/*node*/
	 node* node;

	/*next pointer*/
	struct _list_of_lists* next;

}list_of_lists;



double calc_Q(const int* s, matrix* B, int n);


int division_to_2(matrix* Bg, int* s);

/*double calculate_deltaQ(int* s, spmat* B);*/

node* arry_to_list(int* array, int n);

void list_to_array(node* list, int* array);

void add_group(list_of_lists* groups, node* group);


void unmoved_start(int* unmoved,int ng);

int calc_ng(matrix* B);

void modularity_maximization(matrix* BgHat , int* s);

list_of_lists* divide_network(spmat* A, int size, int* ranks, double* ranks_m);

node* remove_group(list_of_lists* groups);

int is_empty(list_of_lists* groups);

double calc_B_eigen_pair(matrix* B, double *eigenvector, int n);




#endif /* MODULE_ALG_H_ */
