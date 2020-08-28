/*
 * module_alg.c
 *
 *  Created on: 18 באוג׳ 2020
 *      Author: irist
 */

#include "spmat.h"
#include "eigen_pair.h"
#include "B_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>
#include "module_alg.h"


int mult_vectors(const double *v1, const double *v2, int n){
	int i;
	double result = 0.0;
	for(i = 0; i < n; i++){
		result += (*(v1+i))*(*(v2+i));
	}
	return result;
}

double calc_Q(int* s, spmat* B, int n){
	double q;
	double *result = (double*)malloc(n*sizeof(double));
	mult_vector_with_matrix(B, s, result);
	q = 0.5*mult_vectors(result,s);
	free(result);
	return q;
}

int division_to_2(spmat* B, int* g, int n){
	double eigen_val;
	double *eigenvector;
	int *s, i, is_divisible;
	double max_q, q;
	int max_ind;
	is_divisible = 1;
	eigenvector = (double*)malloc(n*sizeof(double));
	s = (int*)malloc(n*sizeof(int));
	spmat *Bg = (spmat*)malloc(sizeof(spmat));
	calc_Bg(B, Bg, g, n);
	calc_Bg_bar(Bg, n);
	eigen_val = calc_B_eigen_pair(Bg, eigenvector, n);
	if(eigen_val <= 0){//The group is indivisible
		is_divisible = 0;
	}
	else{
		for(i = 0; i < n; i++){
			if(eigenvector[i] > 0){
				s[i] = 1;
			}
			else{
				s[i] = -1;
			}
		}
	}
	if(calc_Q(s, B, n) <= 0){//The group is indivisible
		is_divisible = 0;
	}
	else{
		for(i = 0; i < n; i++){
			if(g[i] != 0){
				g[i] = s[i];
			}
		}
	}
	max_q = 0.0;
	max_ind = 0;
	for (i = 0; i < n; i++){
		if(g[i] == -1){
			g[i] = 1;
			q = calc_Q(g, B, n);
			if(q > max_q){
				max_q = q;
				max_ind = i;
				g[i] = -1;
			}
		}
	}
	g[max_ind] = 1;
	free(eigenvector);
	free(s);
	free_in_list(Bg);
	return is_divisible;
}

void add_group(node** groups, int* group, int n){
	int i, j, first;
	node *curr, *elem;
	i = 0;
	first = 1;
	while(groups[i] != NULL){
		i++;
	}
	for(j = 0; j < n; j++){
		if(group[j] != 0){
			elem = (node*)malloc(1*sizeof(node));
			elem->val = group[j];
			elem->next = NULL;
			if(first){
				curr = elem;
				groups[i] = curr;
				first = 0;
			}
			else{
				curr->next = elem;
				curr= curr->next;
				curr->next = NULL;
			}
		}
	}
}

void remove_group(node** groups, int* group){
	int i;
	node *elem;
	i = 0;
	while(groups[i] == NULL){
		i++;
	}
	elem = groups[i];
	while(elem != NULL){
		group[elem->val] = 1;
	}
}

int is_empty(node** groups, int n){
	int i;
	i = 0;
	while((i < n) && (groups[i] == NULL)){
		i++;
	}
	if(i == n){
		return 1;
	}
	else{
		return 0;
	}
}

node** devide_network(spmat* B, int* g, int n){
	int i, j, is_divisible, *g, *g1, *g2;
	int curr_ind_O, curr_ind_P;
	node **groups, **non_divisible_groups, curr;
	groups = (node**)malloc(n*sizeof(node*));
	non_divisible_groups = (node**)malloc(n*sizeof(node*));
	g = (int*)calloc(n, sizeof(int));
	/*create the first group*/
	for(i = 0; i < n; i++){
		g[i] = 1;
	}
	is_divisible = division_to_2(B, g, n);
	if(!is_divisible){
		/*add g to O*/
		add_group(non_divisible_groups, g, n);
	}
	else{
		/*add the 2 groups to P*/
		g1 = (int*)calloc(n, sizeof(int));
		g2 = (int*)calloc(n, sizeof(int));
		for(j = 0; j < n; j++){
			if(g[j] == 1){
				g1[j] = 1;
			}
			else if(g[j] == -1){
				g2[j] = -1;
			}
		}
		add_group(groups, g1, n);
		add_group(groups, g2, n);
		free(g1);
		free(g2);
	}
	free(g);
	while(!is_empty(groups, n)){
		/*remove a group from P and represent the group as an int array g*/
		g = (int*)calloc(n, sizeof(int));
		remove_group(groups, g);
		is_divisible = division_to_2(B, g, n);
			if(!is_divisible){
				/*add g to O*/
				add_group(non_divisible_groups, g, n);
			}
			else{
				/*if one of the groups is empty or has only one element*/




				/*else*/
				/*add the 2 groups to P*/
				g1 = (int*)calloc(n, sizeof(int));
				g2 = (int*)calloc(n, sizeof(int));
					for(j = 0; j < n; j++){
						if(g[j] == 1){
							g1[j] = 1;
						}
						else if(g[j] == -1){
							g2[j] = -1;
						}
					}
					add_group(groups, g1, n);
					add_group(groups, g2, n);
					free(g1);
					free(g2);
			}
			free(g);
	}
	return non_divisible_groups;
}

