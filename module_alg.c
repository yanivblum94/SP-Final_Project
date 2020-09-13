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



double calc_Q(int* s, matrix* B, int n){
	double q;
	double *result = (double*)malloc(n*sizeof(double));
	mult_vector_with_matrix(B, s, result);
	q = 0.5*mult_vectors_int(result, s, n);
	free(result);
	return q;
}

int division_to_2(matrix* Bg, int* s){
	double eigen_val;
	double *eigenvector;
	int size, i, is_divisible;
	size = Bg->size;
	is_divisible = 1;
	eigenvector = (double*)malloc(size*sizeof(double));
	eigen_val = calc_B_eigen_pair(Bg, eigenvector, size);
	if(eigen_val <= 0){/*The group is indivisible*/
		is_divisible = 0;
	}
	else{
		for(i = 0; i < size; i++){
			if(eigenvector[i] > 0){
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

node* arry_to_list(int* array, int n){
	node *curr_node, *elem, *list;
	int i, first;
	first = 1;
	list = (node*)malloc(sizeof(node));
	for(i = 0; i < n; i++){
		if(array[i] != 0){
			elem = (node*)malloc(sizeof(node));
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
	return list;
}

void list_to_array(node* list, int* array){
	while(list != NULL){
		array[list->val] = 1;
		list = list->next;
	}
}


void add_group(list_of_lists* groups, node* group){
	list_of_lists *new_group;
	if(groups->node == NULL){
		new_group = (list_of_lists*)malloc(sizeof(list_of_lists));
		new_group->node = group;
		groups = new_group;
	}
	else{
		while(groups->next != NULL){
				groups = groups->next;
			}
			new_group = (list_of_lists*)malloc(sizeof(list_of_lists));
			new_group->node = group;
			new_group->next = NULL;
			groups->next = new_group;
	}

}

node* remove_group(list_of_lists* groups){
	list_of_lists* elem_to_remove;
	node *group = groups->node;
	elem_to_remove = groups;
	groups = groups->next;
	free(elem_to_remove);
	return group;
}

int is_empty(list_of_lists* groups){
	if(groups == NULL){
		return 1;
	}
	else{
		return 0;
	}
}

list_of_lists* divide_network(spmat* A, int size, int* ranks, double* ranks_m){
	int i, j, is_divisible, *g, *g1, *g2, *s;
	node *group;
	matrix *Bg;
	list_of_lists *groups, *non_divisible_groups;
	char *func = "divide_network";
	groups = (list_of_lists*)calloc(1, sizeof(list_of_lists));
	non_divisible_groups = (list_of_lists*)malloc(sizeof(list_of_lists));
	g = (int*)calloc(size, sizeof(int));
	/*create the first group*/
	for(i = 0; i < size; i++){
		g[i] = 1;
	}
	Bg = allocate_matrix(A, size, ranks, ranks_m, g);
	Bg->c = calc_norm_1(Bg);
	printf("norm is: %f \n" , Bg->c);
	group = arry_to_list(g, size);
	printf("size in matrix  %d \n", Bg->size);
	add_group(groups, group);
	printf("added 1st group\n");
	free(g);
	free(group);
	while(!is_empty(groups)){
		/*remove a group from P and represent the group as an int array g*/
		group = remove_group(groups);
		g = (int*)calloc(size, sizeof(int));
		list_to_array(group, g);
		Bg->g = g;
		s = (int*)calloc(size, sizeof(int));
		is_divisible = division_to_2(Bg, s);
			if(!is_divisible){
				/*add g to O*/
				add_group(non_divisible_groups, group);
				free(group);
			}
			else{
				node *group1, *group2;
				g1 = (int*)calloc(size, sizeof(int));
				g2 = (int*)calloc(size, sizeof(int));
				for(j = 0; j < size; j++){
					if(s[j] == 1){
						g1[j] = 1;
					}
					else if(s[j] == -1){
						g2[j] = -1;
					}
				}
				forceStop(func, 178);
				group1 = arry_to_list(g1, size);
				forceStop(func, 179);
				group2 = arry_to_list(g2, size);
				/*if one of the groups is empty or has only one element*/
				if((group1 == NULL) || (group1->next == NULL)){ /*size of 0 or 1*/
					add_group(non_divisible_groups, group1);/*add to O*/
				}
				else{
					add_group(groups, group1);/*add to P*/
				}
				if((group2 == NULL) || (group2->next == NULL)){ /*size of 0 or 1*/
					add_group(non_divisible_groups, group2);/*add to O*/
				}
				else{
					add_group(groups, group2);/*add to P*/
				}
				free(group1);
				free(group2);
			}
			/*forceStop(func, 195);*/
	}
	free(s);
	forceStop(func, 196);
	return non_divisible_groups;
}

/* calculate Q: by definition
double calculate_deltaQ(int* s, spmat* B){
	double result, *temp=(double*)calloc(1,sizeof(B->n));
	mult_vector_with_matrix(B,s,temp);
	result = mult_vectors(temp, s, B->n);
	return (result*0.5);


}*/

/* initial array values to -1 */
/*void indices_start(int* indices,int n){
	int i;
	for (i = 0; i < n; ++i) {
		indices[i]=-1;
	}
}
/* initial array unmoved to 1 if vertex on g */
void unmoved_start(int* unmoved,int ng){
	int i;
	for (i = 0; i < ng; ++i) {	/*making the unmoved group represnted by array*/
				unmoved[i]=i;
		}
}
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

void modularity_maximization(matrix* B , int* s){
	double *score , *improve ;
	double Q0 , maxscore=0.0, maxImprove, deltaQ;
	int n ,ng, i, j , maxScoreVertex=0, maxImproveIndex;
	int *unmoved, *indices;
	n = B->size;
	deltaQ = calc_Q(s, B, n)*2;
	ng = calc_ng(B);
	unmoved=(int*)calloc(ng,sizeof(int));
	indices=(int*)calloc(ng,sizeof(int));
	score=(double*)calloc(ng,sizeof(double));
	improve=(double*)calloc(ng,sizeof(double));
	while(deltaQ > 0){	/* main while according to line 31 of the alg'*/
	unmoved_start(unmoved, ng);
	for (i = 0; i < ng; ++i) {	/* lines 3-20 alg4*/
			Q0 = calc_Q(s, B, n)*2;
			for (j = 0; j < ng; ++j) {/*for lines 6-10 on alg4 */
					s[j] = -s[j];
					score[j] = (calc_Q(s, B, n)*2)-Q0;
					s[j] = -s[j];
					if(score[j] > maxscore){
						maxScoreVertex = j;
						maxscore = score[j];
					}
			}
			s[maxScoreVertex] = -s[maxScoreVertex];
			indices[i] = maxScoreVertex;
			if(i == 0){
				improve[i] = score[maxScoreVertex];
				maxImprove = improve[i];
				maxImproveIndex = i;
			}else{
				improve[i] = improve[i-1] + score[maxScoreVertex];
				if(improve[i] > maxImprove){ /*maintain maximprove for next part*/
					maxImprove = improve[i];
					maxImproveIndex = i;
				}
			}
			unmoved[maxScoreVertex] = 0;
		}
	for(i = ng-1; i > maxImproveIndex; --i){
		s[indices[i]] = -s[indices[i]];
	}
	if(maxImproveIndex == (ng-1)){
		deltaQ = 0;
	}else{
		deltaQ = improve[maxImproveIndex];
		}
	}
}

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

