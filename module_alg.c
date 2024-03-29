/*
 * module_alg.c
 *
 *      Author: irist
 */

#include "spmat.h"
#include "matrix.h"
#include "eigen_pair.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>
#include "module_alg.h"
#define ISPOSITIVE(X) ((X) > 0.00001)


void free_linked_lists(list_of_lists* sets){
	list_of_lists *temp_set;
	node *list, *temp_list;
	if(sets == NULL){
		free(sets);
	}
	else{
		while(sets != NULL){
			list = sets->node;
			while(list != NULL){
				temp_list = list;
				list = list->next;
				free(temp_list);
			}
			free(list);
			temp_set = sets;
			sets = sets->next;
			free(temp_set);
		}
		free(sets);
	}
}

void free_node(node* group){
	node *temp;
	while(group != NULL){
		temp = group;
		group = group->next;
		free(temp);
	}
	free(group);
}

void print_array(int* g, int n){
	int i;
	for(i = 0; i < n; i++){
		printf("g[%d] = %d\n", i, g[i]);
	}
	printf("done\n");
}
/**
 * terminates the program
 */
void infinite_loop(){
	printf("infinte loop error\n");
	exit(1);
}

/**
 * calculates Q value
 * @param s - s vector
 * @param B - B matrix
 * @param n - size of A
 * @return - Q
 */
double calc_Q(const int* s, matrix* B, int n){
	double q;
	double *result = (double*)calloc(n, sizeof(double));
	mult_matrix_with_vector(B, s, result);
	q = 0.5*(double)(mult_vectors_int(result, s, n));
	free(result);
	return q;
}

/**
 * Algorithm 2 from the doc
 * @param Bg - Bg matrix
 * @param s - s vector
 * @return - whether the group is divisible or not
 */
int division_to_2(matrix* Bg, int* s){
	double eigen_val;
	double *eigenvector;
	int size, i, is_divisible;
	size = Bg->size;
	is_divisible = 1;
	eigenvector = (double*)calloc(size, sizeof(double));
	eigen_val = calc_B_eigen_pair(Bg, eigenvector, size);
	if(eigen_val <= 0){/*The group is indivisible*/
		is_divisible = 0;
	}
	else{
		for(i = 0; i < size; i++){
			if(eigenvector[i] == 0){
				s[i] = 0;
			}
			else if(eigenvector[i] > 0){
				s[i] = 1;
			}
			else{
				s[i] = -1;
			}
		}
	}
	if(calc_Q(s, Bg, size) <= 0){/*The group is indivisible*/
		is_divisible = 0;
	}
	else{
		modularity_maximization(Bg ,s);
	}
	free(eigenvector);
	return is_divisible;
}

/**
 * makes an array into a list
 * @param array - the array we want to make into a list
 * @param n - size of A
 * @return the first node of the list
 */
node* arry_to_list(int* array, int n){
	node *curr_node, *elem, *list;
	int i, first, flag;
	first = 1;
	flag = 1;
	for(i = 0; i < n; i++){
		if(array[i] != 0){
			flag = 0;
			elem = (node*)calloc(1, sizeof(node));
			elem->val = i;
			elem->next = NULL;
			if(first){
				list = elem;
				curr_node = elem;
				first = 0;
			}
			else{
				curr_node->next = elem;
				curr_node = curr_node->next;
				curr_node->next = NULL;
			}
		}
	}
	if(flag){
		list = NULL;
	}
	return list;
}
/**
 * makes a list into an array
 * @param list - the list we get
 * @param array - the array we get in the end
 */
void list_to_array(node* list, int* array){
	while(list != NULL){
		array[list->val] = 1;
		list = list->next;
	}
}


/**
 * adds a group to the division
 * @param groups - the list of groups
 * @param group - the set we want to add
 */
void add_group(list_of_lists* groups, node* group){
	if(groups->node == NULL){/* empty set */
		groups->node = group;
		groups->next = NULL;
	}
	else{
		while(groups->next != NULL){
				groups = groups->next;
			}
		if(groups->next == NULL){
			list_of_lists* new_group;
			new_group = (list_of_lists*)calloc(1, sizeof(list_of_lists));
			new_group->node = group;
			new_group->next = NULL;
			groups->next = new_group;
		}
		else{
			groups->next->node = group;
			groups->next->next = NULL;
		}
	}
}

/**
 * removes a set from the list of sets
 * @param groups - the list of groups
 * @return - reeturns the removed set
 */
node* remove_group(list_of_lists* groups){
	list_of_lists *temp, *elem_to_remove;
	node *group = groups->node;
	if (groups->next == NULL){/* the set contains only one group*/
		groups->node = NULL;
	}
	else{
		temp = groups->next->next;
		elem_to_remove = groups->next;
		groups->node = groups->next->node;
		groups->next = temp;
		free(elem_to_remove);
	}
	return group;
}

/**
 * checks if the list of lists is empty
 * @param groups - the list of sets
 * @return 1 if empty, 0 else
 */
int is_empty(list_of_lists* groups){
	if(groups == NULL || groups->node == NULL){
		return 1;
	}
	else{
		return 0;
	}
}

void del_g(matrix* B){
	int i;
	for(i = 0; i < B->size; i++){
		B->g[i] = 0;
	}
}


/**
 * divides the graph into modularity sets
 * @param A - the matrix graph
 * @param size - size of A
 * @param ranks - ranks vector
 * @param ranks_m - vector of Ki/M
 * @return the list of divided sets
 */
list_of_lists* divide_network(spmat* A, int size, int* ranks, double* ranks_m){
	int i, j, is_divisible,counter=0, *g, *g1, *g2, *s;
	node *group;
	matrix *Bg;
	list_of_lists *groups, *non_divisible_groups;
	groups = (list_of_lists*)calloc(1, sizeof(list_of_lists));
	non_divisible_groups = (list_of_lists*)calloc(1, sizeof(list_of_lists));
	g = (int*)calloc(size, sizeof(int));
	/*create the first group*/
	for(i = 0; i < size; i++){
		g[i] = 1;
	}
	Bg = allocate_matrix(A, size, ranks, ranks_m, g);
	Bg->c = calc_norm_1(Bg);
	group = arry_to_list(g, size);
	add_group(groups, group);
	while(!is_empty(groups)){
		if(counter > 1000 * A->n){
					infinite_loop();
				}
		/*remove a group from P and represent the group as an int array g*/
		group = remove_group(groups);
		del_g(Bg);
		list_to_array(group, Bg->g);
		s = (int*)calloc(size, sizeof(int));
		is_divisible = division_to_2(Bg, s);
			if(!is_divisible){
				/*add g to O*/
				add_group(non_divisible_groups, group);
				free(s);
			}
			else{
				node *group1, *group2;
				free_node(group);
				g1 = (int*)calloc(size, sizeof(int));
				g2 = (int*)calloc(size, sizeof(int));
				for(j = 0; j < size; j++){
					if(s[j] == 1){
						g1[j] = 1;
					}
					else if(s[j] == -1){
						g2[j] = 1;
					}
				}
				free(s);
				group1 = arry_to_list(g1, size);
				group2 = arry_to_list(g2, size);
				free(g1);
				free(g2);
				/*if one of the groups is empty or has only one element*/
				if(group1 == NULL || group2 == NULL){/*size of 0*/
					if(group1 == NULL){
						add_group(non_divisible_groups, group2);/*add to O*/
						free(group1);
					}
					else{
						add_group(non_divisible_groups, group1);/*add to O*/
						free(group2);
					}
				}
				else if(group1->next == NULL){/*size of 1*/
					add_group(non_divisible_groups, group1);/*add to O*/
					add_group(groups, group2);/*add to P*/
				}
				else if(group2->next == NULL){/*size of 1*/
					add_group(non_divisible_groups, group2);/*add to O*/
					add_group(groups, group1);/*add to P*/
				}

				else{
					add_group(groups, group2);/*add to P*/
					add_group(groups, group1);/*add to P*/
				}
	}
			counter++;
	}

	free_matrix(Bg);
	free_linked_lists(groups);
	return non_divisible_groups;
}

/**
 * makes the unmoved array for algo 4
 * @param unmoved - the array we want
 * @param ng - size of Bg
 */
void unmoved_start(int* unmoved,int ng){
	int i;
	for (i = 0; i < ng; ++i) {	/*making the unmoved group represnted by array*/
				unmoved[i]=0;
		}
}
/**
 * calculates the size of Bg
 * @param B - B matrix
 * @return size of Bg
 */
int calc_ng(matrix* B){
	int i, res;
	res =0;
	for(i=0; i < B->size; i++){
		if(B->g[i] != 0){
			res++;
		}
	}
	return res;
}

/**
 * rows 3 - 20 in Algo 4
 * @param B - B matrix
 * @param s - s vector
 * @param ng - size of Bg
 * @param improve - improve indices vector
 * @param unmoved - unmoved vector from the algo
 * @param indices - indices vector from the algo
 */
void find_partition_change(matrix *B, int *s, int ng, double *improve, int *unmoved, int *indices) {
	int i, k, j, index_remove,first;
	double max, q0;
	double *score;
	for (i = 0; i < ng; i++) {
		score = (double*) calloc(ng, sizeof(double));
		q0 = calc_Q(s, B, B->size)*2;
		for (k = 0; k < ng; k++) {
			if(!unmoved[k]){
				s[k] = -s[k];
				score[k] = calc_Q(s, B, B->size)*2 - q0;
				s[k] = -s[k];
			}
		}
		first = 1;
		for (j = 0; j < ng; j++) {
			if(!unmoved[j]){
				if(first){
					max = score[j];
					index_remove = j;
					first=0;
				}else{
					if (score[j] > max) {
						max = score[j];
						index_remove = j;
					}
				}
			}
		}
		s[index_remove] = -s[index_remove];
		indices[i] = index_remove;
		if (i == 0){
			improve[i] = max;
		}
		else{
			improve[i] = improve[i - 1] + max;
		}
		unmoved[index_remove]=1;
		free(score);
	}
}
 /**
  * Algo 4
  * @param B - B matrix
  * @param s - s vector
  */
void modularity_maximization(matrix* B , int* s){
	double *improve ;
	double  max_improve=0, deltaQ;
	int ng, i, j, max_improve_index, counter=0;
	int *unmoved, *indices;
	deltaQ = 1.0;
	ng = calc_ng(B);
	unmoved=(int*)calloc(ng,sizeof(int));
	indices=(int*)calloc(ng,sizeof(int));
	improve=(double*)calloc(ng,sizeof(double));
	while(ISPOSITIVE(deltaQ)){	/* main while according to line 31 of the alg'*/
		if(counter >  1000 * B->size){
					infinite_loop();
				}
		unmoved_start(unmoved, ng);
		find_partition_change(B, s, ng, improve, unmoved, indices);
		max_improve = improve[0];
				max_improve_index = 0;

				for (j = 0; j < ng; j++) {
					if (improve[j] > max_improve) {
						max_improve = improve[j];
						max_improve_index = j;
					}
				}
	for(i = ng-1; i > max_improve_index; --i){
		s[indices[i]] = -s[indices[i]];
	}
	if(max_improve_index == (ng-1)){
		deltaQ = 0;
	}else{
		deltaQ = max_improve;
		}
	counter++;
	}
	free(unmoved);
	free(indices);
	free(improve);
}

/**
 * calculates eigen pair for B
 * @param B - B matrix
 * @param eigenvector - vector of eigens
 * @param n - size of A
 * @return the eigen value
 */
double calc_B_eigen_pair(matrix* B, double *eigenvector, int n){
	double *currvector;
	double eigen_val;
	currvector = (double*)calloc(n, sizeof(double));
	initialize_random_vector(currvector, n);
	calc_eigen(B, currvector, eigenvector, n);
	eigen_val = calc_eigen_val(B, eigenvector, n);
	eigen_val -= B->c;
	free(currvector);
	return eigen_val;

}

