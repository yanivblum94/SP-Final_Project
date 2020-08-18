#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "SPBufferset.h"
#include <math.h>
#include <string.h>

spmat* create_B(spmat* A, int* ranks, int m, int size){
	spmat* B;
	double* row;
	int i,j,col;
	B = spmat_allocate(size);
	for(i=0; i<size; i++){
		row = calloc(size, sizeof(double));
		linked_list curr = A->private[i];
		col = curr.col;
		for(j=0; j< size; j++){
			if(j==col){
				row[j] = 1 - ((ranks[j]*ranks[i])/m);
				curr = curr.next;
			}
			else{
				row[j] = 0 - ((ranks[j]*ranks[i])/m);
			}
		}
		B->add_row(row);
	}
	return B;
}
