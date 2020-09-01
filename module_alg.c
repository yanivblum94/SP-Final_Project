/*
 * module_alg.c
 *
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

void arry_to_list(node* list, int* array, int n){
	node* curr_node, elem;
	int i, first;
	first = 1;
	for(i = 0; i < n; i++){
		if(array[i] != 0){
			elem = (node*)malloc(sizeof(node));
			elem->val = i;
			elem->next = NULL;
			if(first){
				curr_node = elem;
				first = 0;
			}
			else{
				curr_node->next = elem;
				curr_node= curr_node->next;
				curr_node->next = NULL;
			}
		}
	}
}

void list_to_array(node* list, int* array){
	while(list != NULL){
		array[list->val] = 1;
		list = list->next;
	}
}


void add_group(list_of_lists* groups, node* group){
	list_of_lists *new_group;
	while(groups->next != NULL){
		groups = groups->next;
	}
	new_group = (list_of_lists*)malloc(sizeof(list_of_lists));
	new_group->node = group;
	new_group->next = NULL;
	groups->next = new_group;
}

node* remove_group(list_of_lists* groups){
	list_of_lists* elem_to_remove;
	node *group = groups->node;
	elem_to_remove = groups;
	*(groups) = groups->next;
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

list_of_lists* devide_network(spmat* B, int n){
	int i, j, is_divisible, *g, *g1, *g2;
	/*int curr_ind_O, curr_ind_P;*/
	int *g;
	list_of_lists *groups, *non_divisible_groups;
	groups = (list_of_lists*)malloc(sizeof(list_of_lists));
	non_divisible_groups = (list_of_lists*)malloc(sizeof(list_of_lists));
	g = (int*)calloc(n, sizeof(int));
	/*create the first group*/
	for(i = 0; i < n; i++){
		g[i] = 1;
	}
	node* group;
	arry_to_list(group, g, n);
	add_group(groups, group);
	free(g);
	free(group);
	while(!is_empty(groups)){
		/*remove a group from P and represent the group as an int array g*/
		node *group = (node)malloc(sizeof(node));
		remove_group(groups, group);
		g = (int*)calloc(n, sizeof(int));
		list_to_array(group, g);
		is_divisible = division_to_2(B, g, n);
			if(!is_divisible){
				/*add g to O*/
				node* group = (node)malloc(sizeof(node));
				arry_to_list(group, g, n);
				add_group(non_divisible_groups, group);
			}
			else{
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
				node* group1 = (node)malloc(sizeof(node));
				node* group2 = (node)malloc(sizeof(node));
				arry_to_list(group1, g1, n);
				arry_to_list(group2, g2, n);
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
			}
	}
	free(g);
	return non_divisible_groups;
}

/* calculate Q: by definition */
double calculate_deltaQ(int* s, spmat* B){
	double result, *temp=(double*)calloc(1,sizeof(B->n));
	multMatrix(B,s,temp);
	result=multVec(s,temp);
	return (result*0.5);


}

/* initial array values to -1 */
void indices_start(int* indices,int n){
	int i;
	for (i = 0; i < n; ++i) {
		indices[i]=-1;
	}
}
/* initial array unmoved to 1 if vertex on g */
int unmoved_start(int* unmoved,int n,int* s){
	int i,ng=0;
	for (i = 0; i < n; ++i) {	/*making the unmoved group represnted by array*/
			if(s[i]!=0){
				unmoved[i]=1;
				++ng;
			}
		}
	return ng;
}

void modularity_maximization(spmat* BgHat , int* s){
	double* score ,improve ;
	double Q0,maxscore,maxImprove,deltaQ=0;
	int n ,ng, i,j ,maxScoreVertex,maxImproveIndex;;
	int* unmoved,indices;
	n=BgHat->n;
	unmoved=(int*)calloc(1,sizeof(n));
	indices=(int*)calloc(1,sizeof(n));
	score=(double*)calloc(1,sizeof(n));
	improve=(double*)calloc(1,sizeof(n));
	indices_start(indices,n);
	ng=unmoved_start(unmoved,n,s);
	while(deltaQ>=0){	/* main while according to line 31 of the alg'*/
	for (i = 0; i < ng; ++i) {	/* lines 3-20 alg4*/
			Q0 = (calculate_deltaQ(s,BgHat)*2);
			s[0]=-s[0];
			score[0]=(calculate_deltaQ(s,BgHat)*2)-Q0;
			maxscore=score[0];
			maxScoreVertex=0;
			s[0]=-s[0];

			for (j = 1; j < n; ++j) {/*for lines 6-10 on alg4 */
				if(unmoved[j]!=0){
					s[j]=-s[j];
					score[j]=(calculate_deltaQ(s,BgHat)*2)-Q0;
					s[j]=-s[j];
					if(score[j]>maxscore){
						maxScoreVertex=j;
						maxscore=score[j];
					}
				}
			}
			s[maxScoreVertex]=-s[maxScoreVertex];
			indices[i]=maxScoreVertex;
			if(i==0){
				improve[i]=score[maxScoreVertex];
				maxImprove=improve[i];
				maxImproveIndex=i;
			}else{
				improve[i]=improve[i-1]+score[maxScoreVertex];
				if(improve[i]>maxImprove){ /*maintain maximprove for next part*/
					maxImprove=improve[i];
					maxImproveIndex=i;
				}
			}
			unmoved[maxScoreVertex]=0;
		}
	for(i=ng-1;i>maxImproveIndex;--i){
		s[indices[i]]=-s[indices[i]];
	}
	if(maxImproveIndex==(ng-1)){
		deltaQ=0;
	}else{
		deltaQ=improve[maxImproveIndex];
		}
	}
}

