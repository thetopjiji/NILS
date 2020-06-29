/*
 * BCP_c_LS.h
 *
 *  Created on: 6 mai 2019
 *      Author: ren
 */


#include "BCP_struct.h"

#ifndef BCP_C_LS_H_
#define BCP_C_LS_H_

int RepairNF1(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF2(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF3(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF4(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int & retbz);
int RepairNF5(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF6(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairMN1(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int & numI, int & numN);
int RepairMN2(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int & numI, int & numN);
//int RepairBig(int & ite, Struc_Sol *C,  int th);

void DFSReconstruct(Struc_Sol *C, Struc_Sol *Chl);
void CLONEReconstruct(Struc_Sol *Father, Struc_Sol *Mother, Struc_Sol *Chl);
void Modified_DFS(Struc_Sol *C, Struc_Sol *Chl);
void M2_DFS(Struc_Sol *C, Struc_Sol *Chl);


int get_bestCB_pop(int value_alb, int num_pop);
double get_entropy_pop(int value_alb, int num_pop);
double get_averagedistance_pop(int value_alb, int num_pop);
void CopySolution(Struc_Sol *SouSol, Struc_Sol *DesSol);
void EvaSol(Struc_Sol *Sol);
void IniSol(Struc_Sol *Sol);
int get_CycD(int label1, int label2);
void get_newwc(Struc_Sol *C, int *wc, int node1, int node2, int &ncbs);
void get_swap(Struc_Sol *C, int node1, int node2);
void update_cbnodes(Struc_Sol *C, int node1, int node2);
void update_swap(Struc_Sol *C, int node1, int node2,int *wc);
void update_TL(int node1, int node2, int ite);
int Compromise_move(Struc_Sol *C, Struc_Sol* LB, int &ite);
void SI_move(Struc_Sol *C, int mnode, int fnode, int step, int & ite);
void copywc(int *souwc, int *deswc, int l);


int find_midu(Struc_Sol *s, int u);
void newfill_CL(int *oldwc, int *newwc, int *bestwc, int num_wc, int &numI, int &numN, int i, int j);

#endif /* BCP_C_LS_H_ */
