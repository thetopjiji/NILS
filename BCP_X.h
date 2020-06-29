/*
 * BCP_X.h
 *
 *  Created on: 7 mars 2019
 *      Author: ren
 */

#ifndef BCP_X_H_
#define BCP_X_H_
	int judge_one(int *oldwc, int *newwc, int num_wc);
	int judge_or(int *otherwc,int *newwc, int num_wc);
	int judge_part(int *oldwc, int *newwc, int num_wc);
	double judge_fe(int *oldwc, int *newwc, int num_wc);
	double judge_cbs(int oldcbs, int newcbs,int *oldwc, int *newwc,int num_wc);
	void PMX(int *Father, int *Mother, int length, int *c1, int *c2);
	void CX(int *Father, int *Mother, int length, int *c1, int *c2);
	void OX(int *Father, int *Mother, int length, int *c1, int *c2);
	void OX2(int *Father, int *Mother, int length, int *c1, int *c2);
	void OX_modified(int *Father, int *Mother, int length, int *c1, int *c2, int *Fcbnodes, int *Mcbnodes, int *PermuF, int *PermuM, int cb);
	int DPX(int *Father, int *Mother, int length, int *c1, int *c2, int **mat_s);
	int get_dis_tsp(int *Father, int *Mother, int length, int **matrix_sum);
#endif /* BCP_X_H_ */
