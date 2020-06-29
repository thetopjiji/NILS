/*
 * BCP_struct.h
 *
 *  Created on: 6 mai 2019
 *      Author: ren
 */

#ifndef BCP_STRUCT_H_
#define BCP_STRUCT_H_


typedef struct struct_individual {
	int *permutation;
	int *permutationNew;
	int cbmp;
	int cbs;
	int *wc;
	int *cbnodes;
} Struc_Sol;
//Struc_Sol *CurrentS,*BestS, *PopS, *ChS1, *ChS2;




#endif /* BCP_STRUCT_H_ */
