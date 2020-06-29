/*
 * BCP_c_LS.cpp
 *
 *  Created on: 6 mai 2019
 *      Author: ren
 */


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <cmath>
#include <math.h>
#include <complex>
#include <assert.h>
#include <iomanip>
#include <sys/times.h>
#include <unistd.h>
#include "BCP_struct.h"
#include "BCP_X.h"

using namespace std;
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

//extern variable
extern int num_vex;
//extern int num_unk;
extern int num_edge;
//extern int max_deg;
extern int **start_end;
extern int **edge;
//extern int *degree_node;
extern int **simple_edge;
extern int alreadybest;

//struct tms glo_start;
//struct tms glo_end;
//struct tms midt;
//double clockTicksPerSecond;
//double startTimeSeconds;
//double cpuTimeTheBest;

extern Struc_Sol *CurrentS, *PopS;
//extern Struc_Sol *BestS,*ChS1,*ChS2,*LBS,*lb_local,*llb,*petlb;
extern int *CV;
extern int **Candidatemove;
extern int **mat_sum;
extern int **TL;
extern int **CandiList;
extern int **CandiSI;

//extern int len_ls;
//extern double len_th;
//extern double p_best;
//extern double p_NF;
extern double percent_can;
//extern int depth;
//extern int max_cbs;

extern double alpha;

//local variable
int list_aj[15]={1,2,1,4,1,2,1,8,1,2,1,4,1,2,1};

//declaration of function

int RepairNF1(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF2(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF3(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF4(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int &retbz);
int RepairNF5(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);
int RepairNF6(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);

int RepairMN1(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int & numI, int & numN);
int RepairMN2(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int & numI, int & numN);

void DFSReconstruct(Struc_Sol *C, Struc_Sol *Chl);
void CLONEReconstruct(Struc_Sol *Father, Struc_Sol *Mother, Struc_Sol *Chl);
void Modified_DFS(Struc_Sol *C, Struc_Sol *Chl);
void M2_DFS(Struc_Sol *C, Struc_Sol *Chl);


//int RepairBig(int & ite, Struc_Sol *C, Struc_Sol *LO, int th);

int shift_insert(Struc_Sol *C);
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
int calcul_tenure(int ite);
int Compromise_move(Struc_Sol *C, Struc_Sol* LB, int &ite);
void SI_move(Struc_Sol *C, int mnode, int fnode, int step, int & ite);
int fix_label(int inputlabel);
void copywc(int *souwc, int *deswc, int l);
//local declaration
int find_midu(Struc_Sol *s, int u);

void fill_CL(int *oldwc, int *newwc, int *bestwc, int num_wc, int &nimp, int i, int j);
void newfill_CL(int *oldwc, int *newwc, int *bestwc, int num_wc, int &numI, int &numN, int i, int j);

void generate_sets(int *set_s, int l);
int find_balance(Struc_Sol *C,int *wc,int node1, int node2);
void balance_CL(int i, int j);
void cbs_Candidate(int *oldwc, int *newwc, int *bestwc, int num_wc,int oldcbs, int newcbs, int &bestcbs,  int &nimp, int i, int j);

int choose_one_suit(int *ArrCV, int *ArrRT, int flag);
void arr_node(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew);
void arr_fillNext(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr);
int get_realpos(int fixlabel, int step);

void arr_node_clone(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, Struc_Sol *M);
void arr_fillNext_clone(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr,Struc_Sol *M);

void Modarr_node(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew);
void Modarr_fillNext(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr);

void M2arr_node(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew);
void M2arr_fillNext(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr);


//function
void M2arr_fillNext(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr){
	int i,j;
	int stnode, ennode;
//	int interns, interne;
	int tempnode;
	int startpos1;
	int startpos2;
	int num_can;
	int ssq[num_vex];
	int Nchoi;

	int realpos;
	int fixlabel=PermC[ro];

	ArrRT[ro]=1;

	stnode=start_end[ro][0];
	ennode=start_end[ro][1];

	for(i=stnode;i<=ennode;i++){
		tempnode=edge[i][1];
		num_can=0;
		Nchoi=-1;
		if(ArrCV[tempnode]==0 && ArrRT[tempnode]==1){
			cout<<"error in Modarr_fillNext, not set but as root"<<endl;
			exit(-3);
		}
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==1) continue;		//already set and as root
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==0) {				//already set but not as root
			Next_root[num_nr]=tempnode;
			num_nr++;
			ArrRT[tempnode]=1;

		}
		else if(ArrCV[tempnode]==0 && ArrRT[tempnode]==0){				//not set and not as root
			startpos1=get_realpos(fixlabel,alreadybest);
			startpos2=get_realpos(fixlabel,-alreadybest);

			for(j=0;j<=alreadybest;j++){				//arrange the pos
				realpos=get_realpos(startpos1,-j);
				if(PermCnew[realpos]==-1){
					ssq[num_can]=realpos;
					num_can++;
				}
				realpos=get_realpos(startpos2,j);
				if(PermCnew[realpos]==-1){
					ssq[num_can]=realpos;
					num_can++;
				}
			}

			if(num_can>0){
				Nchoi=rand()%num_can;
				realpos=ssq[Nchoi];
				PermCnew[realpos]=tempnode;		//update the PermCnew
				PermC[tempnode]=realpos;		//update the PermC
				ArrCV[tempnode]=1;				//update the CV

			}
			else{
				for(j=0;j<=num_vex;j++){				//arrange the pos
					realpos=get_realpos(startpos1,j);
					if(PermCnew[realpos]==-1){
						PermCnew[realpos]=tempnode;		//update the PermCnew
						PermC[tempnode]=realpos;		//update the PermC
						ArrCV[tempnode]=1;				//update the CV
						break;
					}
					realpos=get_realpos(startpos2,-j);
					if(PermCnew[realpos]==-1){
						PermCnew[realpos]=tempnode;
						PermC[tempnode]=realpos;
						ArrCV[tempnode]=1;
						break;
					}
				}
			}
			if(j>=num_vex+1){
				cout<<"error in Modarr_fillNext, no space"<<endl;
				exit(-1);
			}
			Next_root[num_nr]=tempnode;	//update the Next_root
			num_nr++;
			ArrRT[tempnode]=1;

		}
		else {
			cout<<"error in arr_fillNext, other condition"<<endl;
			exit(-3);
		}

	}
	if(num_nr>=num_vex){
	  cout<<"error in num_nr"<<endl;
	  exit(-1);
	}

}

void M2arr_node(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew){
	int i;
	int Cur_root[num_vex];
	int Next_root[num_vex];
	int num_cr=-1;
	int num_nr=-1;

	/*judge whether the rootN is set*/
	if(ArrCV[rootN]==0) {
		for(i=0;i<num_vex;i++){
			if(PermCnew[i]==-1){
				PermCnew[i]=rootN;
				PermC[rootN]=i;
				CV[rootN]=1;
				break;
			}
		}
		if(i==num_vex){
			cout<<"error in arraging label to unset node"<<endl;
			exit(-1);
		}
	}
	Next_root[0]=rootN;
	num_nr=1;
	ArrRT[rootN]=1;					// rootN is written in Next_root
	/*start to explore the tree*/
	do{
		copywc(Next_root,Cur_root,num_nr);
		num_cr=num_nr;
		num_nr=0;
		/*From Cur_root to explore*/
		for(i=0;i<num_cr;i++)
			M2arr_fillNext(Cur_root[i], ArrCV, ArrRT, PermC, PermCnew, Cur_root, Next_root, num_nr);
	}while(num_nr!=0);
	return;
}

void M2_DFS(Struc_Sol *C, Struc_Sol *Chl){
	int i;
	int cb=C->cbmp;
	int Asroot[num_vex];						// already as root
	int rootN=-1;
	int PermC[num_vex];
	int PermCnew[num_vex];
	int flag;

	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]<=alpha*cb) CV[i]=0;		// will change
		else CV[i]=1;							// dont change
		Asroot[i]=0;
		PermC[i]=-1;
		PermCnew[i]=-1;
	}
	/*The initialisation of Childsolution*/
	for(i=0;i<num_vex;i++){
		if(CV[i]==1){
			PermC[i]=C->permutation[i];
			PermCnew[PermC[i]]=i;
		}
	}

	/************Get the rootN only set************/
	rootN=choose_one_suit(CV,Asroot,1);
	if(rootN==-1){
		rootN=choose_one_suit(CV,Asroot,0);
		M2arr_node(rootN,CV,Asroot,PermC,PermCnew);					//all the nodes get to reset
	}
	else M2arr_node(rootN,CV,Asroot,PermC,PermCnew);
	/************For the case of there are nodes set but not as root****************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]+Asroot[i]==1) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,2);
		M2arr_node(rootN,CV,Asroot,PermC,PermCnew);
	}while(flag);
	/**********For the case of the nodes which are not set**************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]==0) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,0);
		M2arr_node(rootN,CV,Asroot,PermC,PermCnew);
	}while(flag);

//	/*check the result*/
//	int check_seq[num_vex];
//	for(i=0;i<num_vex;i++) check_seq[i]=0;
//	for(i=0;i<num_vex;i++){
//		if(check_seq[PermC[i]]) {
//			cout<<"error in DFS PermC, duplicate"<<endl;
//			exit(-2);
//		}
//		check_seq[PermC[i]]=1;
//	}
//
//	for(i=0;i<num_vex;i++) check_seq[i]=0;
//	for(i=0;i<num_vex;i++){
//		if(check_seq[PermCnew[i]]) {
//			cout<<"error in DFS PermCnew, duplicate"<<endl;
//			exit(-2);
//		}
//		check_seq[PermCnew[i]]=1;
//	}
//
//	for(i=0;i<num_vex;i++){
//		if(PermC[PermCnew[i]]!=i){
//			cout<<"error in DFS, not according"<<endl;
//			exit(-2);
//		}
//	}

	/*Copy and return*/
	copywc(PermC,Chl->permutation,num_vex);
	copywc(PermCnew,Chl->permutationNew,num_vex);

}



void Modarr_fillNext(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr){
	int i,j;
	int stnode, ennode;
//	int interns, interne;
	int tempnode;
	int startpos1;
	int startpos2;

	int realpos;
	int fixlabel=PermC[ro];

	ArrRT[ro]=1;

	stnode=start_end[ro][0];
	ennode=start_end[ro][1];

	for(i=stnode;i<=ennode;i++){
		tempnode=edge[i][1];
		if(ArrCV[tempnode]==0 && ArrRT[tempnode]==1){
			cout<<"error in Modarr_fillNext, not set but as root"<<endl;
			exit(-3);
		}
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==1) continue;		//already set and as root
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==0) {				//already set but not as root
			Next_root[num_nr]=tempnode;
			num_nr++;
			ArrRT[tempnode]=1;

		}
		else if(ArrCV[tempnode]==0 && ArrRT[tempnode]==0){				//not set and not as root
			startpos1=get_realpos(fixlabel,alreadybest);
			startpos2=get_realpos(fixlabel,-alreadybest);

			for(j=0;j<=num_vex;j++){				//arrange the pos
				realpos=get_realpos(startpos1,-j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;		//update the PermCnew
					PermC[tempnode]=realpos;		//update the PermC
					ArrCV[tempnode]=1;				//update the CV
					break;
				}
				realpos=get_realpos(startpos2,j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;
					PermC[tempnode]=realpos;
					ArrCV[tempnode]=1;
					break;
				}
			}
			if(j>=num_vex+1){
				cout<<"error in Modarr_fillNext, no space"<<endl;
				exit(-1);
			}
			Next_root[num_nr]=tempnode;	//update the Next_root
			num_nr++;
			ArrRT[tempnode]=1;

		}
		else {
			cout<<"error in arr_fillNext, other condition"<<endl;
			exit(-3);
		}

	}
	if(num_nr>=num_vex){
	  cout<<"error in num_nr"<<endl;
	  exit(-1);
	}

}

void Modarr_node(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew){
	int i;
	int Cur_root[num_vex];
	int Next_root[num_vex];
	int num_cr=-1;
	int num_nr=-1;

	/*judge whether the rootN is set*/
	if(ArrCV[rootN]==0) {
		for(i=0;i<num_vex;i++){
			if(PermCnew[i]==-1){
				PermCnew[i]=rootN;
				PermC[rootN]=i;
				CV[rootN]=1;
				break;
			}
		}
		if(i==num_vex){
			cout<<"error in arraging label to unset node"<<endl;
			exit(-1);
		}
	}
	Next_root[0]=rootN;
	num_nr=1;
	ArrRT[rootN]=1;					// rootN is written in Next_root
	/*start to explore the tree*/
	do{
		copywc(Next_root,Cur_root,num_nr);
		num_cr=num_nr;
		num_nr=0;
		/*From Cur_root to explore*/
		for(i=0;i<num_cr;i++)
			Modarr_fillNext(Cur_root[i], ArrCV, ArrRT, PermC, PermCnew, Cur_root, Next_root, num_nr);
	}while(num_nr!=0);
	return;
}

void Modified_DFS(Struc_Sol *C, Struc_Sol *Chl){
	int i;
	int cb=C->cbmp;
	int Asroot[num_vex];						// already as root
	int rootN=-1;
	int PermC[num_vex];
	int PermCnew[num_vex];
	int flag;

	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]<=0.8*cb) CV[i]=0;		// will change
		else CV[i]=1;							// dont change
		Asroot[i]=0;
		PermC[i]=-1;
		PermCnew[i]=-1;
	}
	/*The initialisation of Childsolution*/
	for(i=0;i<num_vex;i++){
		if(CV[i]==1){
			PermC[i]=C->permutation[i];
			PermCnew[PermC[i]]=i;
		}
	}

	/************Get the rootN only set************/
	rootN=choose_one_suit(CV,Asroot,1);
	if(rootN==-1){
		rootN=choose_one_suit(CV,Asroot,0);
		Modarr_node(rootN,CV,Asroot,PermC,PermCnew);					//all the nodes get to reset
	}
	else Modarr_node(rootN,CV,Asroot,PermC,PermCnew);
	/************For the case of there are nodes set but not as root****************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]+Asroot[i]==1) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,2);
		Modarr_node(rootN,CV,Asroot,PermC,PermCnew);
	}while(flag);
	/**********For the case of the nodes which are not set**************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]==0) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,0);
		Modarr_node(rootN,CV,Asroot,PermC,PermCnew);
	}while(flag);

//	/*check the result*/
//	int check_seq[num_vex];
//	for(i=0;i<num_vex;i++) check_seq[i]=0;
//	for(i=0;i<num_vex;i++){
//		if(check_seq[PermC[i]]) {
//			cout<<"error in DFS PermC, duplicate"<<endl;
//			exit(-2);
//		}
//		check_seq[PermC[i]]=1;
//	}
//
//	for(i=0;i<num_vex;i++) check_seq[i]=0;
//	for(i=0;i<num_vex;i++){
//		if(check_seq[PermCnew[i]]) {
//			cout<<"error in DFS PermCnew, duplicate"<<endl;
//			exit(-2);
//		}
//		check_seq[PermCnew[i]]=1;
//	}
//
//	for(i=0;i<num_vex;i++){
//		if(PermC[PermCnew[i]]!=i){
//			cout<<"error in DFS, not according"<<endl;
//			exit(-2);
//		}
//	}

	/*Copy and return*/
	copywc(PermC,Chl->permutation,num_vex);
	copywc(PermCnew,Chl->permutationNew,num_vex);

}




void arr_fillNext_clone(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr,Struc_Sol *M){
	int i,j;
	int stnode, ennode;
//	int interns, interne;
	int tempnode;

	int realpos;
	int fixlabel=PermC[ro];
	int ano_len;
	int startpos1;
	int startpos2;

	ArrRT[ro]=1;

	stnode=start_end[ro][0];
	ennode=start_end[ro][1];

	for(i=stnode;i<=ennode;i++){
		tempnode=edge[i][1];

		if(ArrCV[tempnode]==0 && ArrRT[tempnode]==1){
			cout<<"error in arr_fillNext, not set but as root"<<endl;
			exit(-3);
		}
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==1) continue;		//already set and as root
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==0) {				//already set but not as root
			Next_root[num_nr]=tempnode;
			num_nr++;
			ArrRT[tempnode]=1;

		}
		else if(ArrCV[tempnode]==0 && ArrRT[tempnode]==0){				//not set and not as root
			ano_len=get_CycD(M->permutation[ro],M->permutation[tempnode]);
			startpos1=get_realpos(fixlabel,ano_len);
			startpos2=get_realpos(fixlabel,-ano_len);

			for(j=0;j<=num_vex/2;j++){				//arrange the pos
				realpos=get_realpos(startpos1,j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;		//update the PermCnew
					PermC[tempnode]=realpos;		//update the PermC
					ArrCV[tempnode]=1;				//update the CV
					break;
				}
				realpos=get_realpos(startpos1,-j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;
					PermC[tempnode]=realpos;
					ArrCV[tempnode]=1;
					break;
				}

				realpos=get_realpos(startpos2,j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;		//update the PermCnew
					PermC[tempnode]=realpos;		//update the PermC
					ArrCV[tempnode]=1;				//update the CV
					break;
				}
				realpos=get_realpos(startpos2,-j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;
					PermC[tempnode]=realpos;
					ArrCV[tempnode]=1;
					break;
				}

			}
			if(j>num_vex/2){
				cout<<"error in arr_fillNext_clone, no rest label"<<endl;
				exit(-3);
			}
			Next_root[num_nr]=tempnode;	//update the Next_root
			num_nr++;
			ArrRT[tempnode]=1;

		}
		else {
			cout<<"error in arr_fillNext_clone, other condition"<<endl;
			exit(-3);
		}

	}
	if(num_nr>=num_vex){
	  cout<<"error in num_nr"<<endl;
	  exit(-1);
	}

}


void arr_node_clone(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, Struc_Sol *M){
	int i;
	int Cur_root[num_vex];
	int Next_root[num_vex];
	int num_cr=-1;
	int num_nr=-1;

	/*judge whether the rootN is set*/
	if(ArrCV[rootN]==0) {
		for(i=0;i<num_vex;i++){
			if(PermCnew[i]==-1){
				PermCnew[i]=rootN;
				PermC[rootN]=i;
				CV[rootN]=1;
				break;
			}
		}
		if(i==num_vex){
			cout<<"error in arraging label to unset node"<<endl;
			exit(-1);
		}
	}
	Next_root[0]=rootN;
	num_nr=1;
	ArrRT[rootN]=1;					// rootN is written in Next_root
	/*start to explore the tree*/
	do{
		copywc(Next_root,Cur_root,num_nr);
		num_cr=num_nr;
		num_nr=0;
		/*From Cur_root to explore*/
		for(i=0;i<num_cr;i++)
			arr_fillNext_clone(Cur_root[i], ArrCV, ArrRT, PermC, PermCnew, Cur_root, Next_root, num_nr,M);
	}while(num_nr!=0);
	return;
}


void CLONEReconstruct(Struc_Sol *Father, Struc_Sol *Mother, Struc_Sol *Chl){
	//make the operation on Solution Father

	int i;
	int cb=Father->cbmp;
	int Asroot[num_vex];						// already as root
	int rootN=-1;
	int PermC[num_vex];
	int PermCnew[num_vex];
	int flag;

	for(i=0;i<num_vex;i++){
		if(Father->cbnodes[i]<=0.8*cb) CV[i]=0;		// will change
		else CV[i]=1;							// dont change
		Asroot[i]=0;
		PermC[i]=-1;
		PermCnew[i]=-1;
	}

	/*The initialisation of Childsolution*/
	for(i=0;i<num_vex;i++){
		if(CV[i]==1){
			PermC[i]=Father->permutation[i];
			PermCnew[PermC[i]]=i;
		}
	}
	/************Get the rootN only set************/
	rootN=choose_one_suit(CV,Asroot,1);
	if(rootN==-1){
		rootN=choose_one_suit(CV,Asroot,0);
		arr_node_clone(rootN,CV,Asroot,PermC,PermCnew,Mother);					//all the nodes get to reset
	}
	else arr_node_clone(rootN,CV,Asroot,PermC,PermCnew,Mother);
	/************For the case of there are nodes set but not as root****************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]+Asroot[i]==1) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,2);
		arr_node_clone(rootN,CV,Asroot,PermC,PermCnew,Mother);
	}while(flag);
	/**********For the case of the nodes which are not set**************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]==0) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,0);
		arr_node_clone(rootN,CV,Asroot,PermC,PermCnew,Mother);
	}while(flag);

//	/*check the result*/
	int check_seq[num_vex];
	for(i=0;i<num_vex;i++) check_seq[i]=0;
	for(i=0;i<num_vex;i++){
		if(check_seq[PermC[i]]) {
			cout<<"error in DFS PermC, duplicate"<<endl;
			exit(-2);
		}
		check_seq[PermC[i]]=1;
	}

	for(i=0;i<num_vex;i++) check_seq[i]=0;
	for(i=0;i<num_vex;i++){
		if(check_seq[PermCnew[i]]) {
			cout<<"error in DFS PermCnew, duplicate"<<endl;
			exit(-2);
		}
		check_seq[PermCnew[i]]=1;
	}

	for(i=0;i<num_vex;i++){
		if(PermC[PermCnew[i]]!=i){
			cout<<"error in DFS, not according"<<endl;
			exit(-2);
		}
	}

	//CopySolution
	copywc(PermC,Chl->permutation,num_vex);
	copywc(PermCnew,Chl->permutationNew,num_vex);

}


int get_realpos(int fixlabel, int step){
	int realpos;
	realpos=(fixlabel+num_vex+step)%num_vex;
	return realpos;
}

void arr_fillNext(int ro, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew, int *Cur_root, int *Next_root, int &num_nr){
	int i,j;
	int stnode, ennode;
//	int interns, interne;
	int tempnode;

	int realpos;
	int fixlabel=PermC[ro];

	ArrRT[ro]=1;

	stnode=start_end[ro][0];
	ennode=start_end[ro][1];

	for(i=stnode;i<=ennode;i++){
		tempnode=edge[i][1];
		if(ArrCV[tempnode]==0 && ArrRT[tempnode]==1){
			cout<<"error in arr_fillNext, not set but as root"<<endl;
			exit(-3);
		}
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==1) continue;		//already set and as root
		else if(ArrCV[tempnode]==1 && ArrRT[tempnode]==0) {				//already set but not as root
			Next_root[num_nr]=tempnode;
			num_nr++;
			ArrRT[tempnode]=1;
//			interns=start_end[tempnode][0];
//			interne=start_end[tempnode][1];
//			for(j=interns;j<=interne;j++){
//				internnode=edge[j][1];
//				if(ArrRT[internnode]==1) continue;
//				Next_root[num_nr]=internnode;		//update Next_root
//				num_nr++;
//			}

		}
		else if(ArrCV[tempnode]==0 && ArrRT[tempnode]==0){				//not set and not as root

			for(j=1;j<=num_vex/2;j++){				//arrange the pos
				realpos=get_realpos(fixlabel,j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;		//update the PermCnew
					PermC[tempnode]=realpos;		//update the PermC
					ArrCV[tempnode]=1;				//update the CV
					break;
				}
				realpos=get_realpos(fixlabel,-j);
				if(PermCnew[realpos]==-1){
					PermCnew[realpos]=tempnode;
					PermC[tempnode]=realpos;
					ArrCV[tempnode]=1;
					break;
				}
			}
			Next_root[num_nr]=tempnode;	//update the Next_root
			num_nr++;
			ArrRT[tempnode]=1;
//			interns=start_end[tempnode][0];
//			interne=start_end[tempnode][1];
//			for(j=interns;j<=interne;j++){
//				internnode=edge[j][1];
//				if(ArrRT[internnode]==1) continue;
//				Next_root[num_nr]=internnode;	//update the Next_root
//				num_nr++;
//			}

		}
		else {
			cout<<"error in arr_fillNext, other condition"<<endl;
			exit(-3);
		}

	}
	if(num_nr>=num_vex){
	  cout<<"error in num_nr"<<endl;
	  exit(-1);
	}

}

void arr_node(int rootN, int *ArrCV, int *ArrRT, int *PermC, int *PermCnew){

	int i;
	int Cur_root[num_vex];
	int Next_root[num_vex];
	int num_cr=-1;
	int num_nr=-1;

	/*judge whether the rootN is set*/
	if(ArrCV[rootN]==0) {
		for(i=0;i<num_vex;i++){
			if(PermCnew[i]==-1){
				PermCnew[i]=rootN;
				PermC[rootN]=i;
				CV[rootN]=1;
				break;
			}
		}
		if(i==num_vex){
			cout<<"error in arraging label to unset node"<<endl;
			exit(-1);
		}
	}
	Next_root[0]=rootN;
	num_nr=1;
	ArrRT[rootN]=1;					// rootN is written in Next_root
	/*start to explore the tree*/
	do{
		copywc(Next_root,Cur_root,num_nr);
		num_cr=num_nr;
		num_nr=0;
		/*From Cur_root to explore*/
		for(i=0;i<num_cr;i++)
			arr_fillNext(Cur_root[i], ArrCV, ArrRT, PermC, PermCnew, Cur_root, Next_root, num_nr);
	}while(num_nr!=0);
	return;
}

int choose_one_suit(int *ArrCV, int *ArrRT, int flag){
	int node=-1;
	int i;
	int tempseq[num_vex];
	int count=0;

	if(flag==1){					// to choose one alreay set
		for(i=0;i<num_vex;i++)
			if(ArrCV[i]==1){
				tempseq[count]=i;
				count++;
			}
	}
	else if(flag==0){				// to choose one not set
		for(i=0;i<num_vex;i++)
			if(ArrCV[i]==0){
				tempseq[count]=i;
				count++;
			}
	}
	else{
		for(i=0;i<num_vex;i++)		// to choose one set but not as root
			if(ArrCV[i]+ArrRT[i]==1){
				tempseq[count]=i;
				count++;
			}
	}

	if(count==0) return node;
	node=tempseq[rand()%count];
	return node;
}

void DFSReconstruct(Struc_Sol *C, Struc_Sol *Chl){
	int i;
	int cb=C->cbmp;
	int Asroot[num_vex];						// already as root
	int rootN=-1;
	int PermC[num_vex];
	int PermCnew[num_vex];
	int flag;

	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]<=0.8*cb) CV[i]=0;		// will change
		else CV[i]=1;							// dont change
		Asroot[i]=0;
		PermC[i]=-1;
		PermCnew[i]=-1;
	}
	/*The initialisation of Childsolution*/
	for(i=0;i<num_vex;i++){
		if(CV[i]==1){
			PermC[i]=C->permutation[i];
			PermCnew[PermC[i]]=i;
		}
	}

	/************Get the rootN only set************/
	rootN=choose_one_suit(CV,Asroot,1);
	if(rootN==-1){
		rootN=choose_one_suit(CV,Asroot,0);
		arr_node(rootN,CV,Asroot,PermC,PermCnew);					//all the nodes get to reset
	}
	else arr_node(rootN,CV,Asroot,PermC,PermCnew);
	/************For the case of there are nodes set but not as root****************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]+Asroot[i]==1) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,2);
		arr_node(rootN,CV,Asroot,PermC,PermCnew);
	}while(flag);
	/**********For the case of the nodes which are not set**************/
	do{
		flag=0;
		for(i=0;i<num_vex;i++) {
			if(CV[i]==0) {
				flag=1;
				break;
			}
		}
		if(flag==0) break;
		rootN=choose_one_suit(CV,Asroot,0);
		arr_node(rootN,CV,Asroot,PermC,PermCnew);
	}while(flag);

//	/*check the result*/
//	int check_seq[num_vex];
//	for(i=0;i<num_vex;i++) check_seq[i]=0;
//	for(i=0;i<num_vex;i++){
//		if(check_seq[PermC[i]]) {
//			cout<<"error in DFS PermC, duplicate"<<endl;
//			exit(-2);
//		}
//		check_seq[PermC[i]]=1;
//	}
//
//	for(i=0;i<num_vex;i++) check_seq[i]=0;
//	for(i=0;i<num_vex;i++){
//		if(check_seq[PermCnew[i]]) {
//			cout<<"error in DFS PermCnew, duplicate"<<endl;
//			exit(-2);
//		}
//		check_seq[PermCnew[i]]=1;
//	}
//
//	for(i=0;i<num_vex;i++){
//		if(PermC[PermCnew[i]]!=i){
//			cout<<"error in DFS, not according"<<endl;
//			exit(-2);
//		}
//	}

	/*Copy and return*/
	copywc(PermC,Chl->permutation,num_vex);
	copywc(PermCnew,Chl->permutationNew,num_vex);

}



int fix_label(int inputlabel){
	return (inputlabel+num_vex)%num_vex;
}

void SI_move(Struc_Sol *C, int mnode, int fnode, int step, int & ite){
	int start_label;
	int stop_label;
	int dir_move=0;

	int nowlabel;
	int nextlabel;
	int nextnode;
	int i;
	/*Decide the dirtion of rotation*/
	start_label=C->permutation[mnode];
	stop_label=C->permutation[fnode];
	if((stop_label>start_label)&&(stop_label-start_label<num_vex/2)) dir_move=1;
	else if((stop_label>start_label)&&(stop_label-start_label>=num_vex/2)) dir_move=-1;
	else if((stop_label<start_label)&&(start_label-stop_label<num_vex/2)) dir_move=-1;
	else if((stop_label<start_label)&&(start_label-stop_label>=num_vex/2)) dir_move=1;
	/*get move*/
	nowlabel=start_label;
	for(i=0;i<step;i++){
		nextlabel=fix_label(nowlabel+dir_move);
		nextnode=C->permutationNew[nextlabel];
		/*adding part*/
//		if(TL[mnode][nextnode]>=ite) break;
//		if(TL[nextnode][mnode]>=ite) break;
		/************************************/
		get_swap(C,mnode,nextnode);
		nowlabel=nextlabel;

		/*update TL*/
//		update_TL(mnode,nextnode,ite);
	}

}


int Compromise_move(Struc_Sol *C, Struc_Sol* LB, int &ite){
	int i;
	int nn_cv=0;
	int tC[num_vex];
	int index_chosen;
	int node_chosen;
	int temp_node;
	int pair_node=-1;
	int ss,ee;
	int largestcb=0;
	int temp_cb;
	int step=0;


	for(i=0;i<num_vex;i++)
		if(C->cbnodes[i]>=C->cbmp*1.0){
			tC[nn_cv]=i;
			nn_cv++;
		}
	index_chosen=rand()%nn_cv;
	node_chosen=tC[index_chosen];
	ss=start_end[node_chosen][0];
	ee=start_end[node_chosen][1];

	for(i=ss;i<=ee;i++){
		temp_node=edge[i][1];
		temp_cb=get_CycD(C->permutation[node_chosen],C->permutation[temp_node]);
		if(temp_cb>largestcb) {
			pair_node=temp_node;
			largestcb=temp_cb;
		}
	}

	if(pair_node==-1) {
		cout<<"error in pairnode"<<endl;
		exit(-1);
	}

	step=rand()%largestcb;
	if(rand()/(RAND_MAX+1.0)>=0.5) SI_move(C,node_chosen,pair_node,step, ite);
	else SI_move(C,pair_node,node_chosen,step,ite);

	/*adding part*/
//	ite++;
	/*****************************/
	return 0;
}


void cbs_Candidate(int *oldwc, int *newwc, int *bestwc, int num_wc,int oldcbs, int newcbs, int &bestcbs,  int &nimp, int i, int j){
	double judge_f;
	judge_f=judge_cbs(oldcbs,newcbs,oldwc,newwc,num_wc);
	if(judge_f<=0){
		if(nimp==0){
			CandiList[0][0]=i;
			CandiList[0][1]=j;
			copywc(newwc,bestwc,num_wc);
			bestcbs=newcbs;
			nimp++;
		}
		else {
			if(judge_cbs(bestcbs,newcbs,bestwc,newwc,num_wc)<0){
//			if(judge_part(bestwc,newwc,num_wc)<0){
				CandiList[nimp][0]=CandiList[0][0];
				CandiList[nimp][1]=CandiList[0][1];
				CandiList[0][0]=i;
				CandiList[0][1]=j;
				copywc(newwc,bestwc,num_wc);
				bestcbs=newcbs;
				nimp++;
			}
			else {
				CandiList[nimp][0]=i;
				CandiList[nimp][1]=j;
				nimp++;
			}
		}
	}
}


void generate_sets(int *set_s, int l){
	int i=0,tt;
	int temp[l];
	for(i=0;i<l;i++) {
		temp[i]=rand()%l;
		set_s[i]=i;
	}
	for(i=0;i<l;i++){
		tt=set_s[i];
		set_s[i]=set_s[temp[i]];
		set_s[temp[i]]=tt;
	}
}

void fill_CL(int *oldwc, int *newwc, int *bestwc, int num_wc, int &nimp, int i, int j){
//	double judge_f;
	int judge_f;
	judge_f=judge_part(oldwc,newwc,num_wc);
	if(judge_f<=0){
		if(nimp==0){
			CandiList[0][0]=i;
			CandiList[0][1]=j;
			copywc(newwc,bestwc,num_wc);
			nimp++;
		}
		else {
			if(judge_part(bestwc,newwc,num_wc)<0){
				CandiList[nimp][0]=CandiList[0][0];
				CandiList[nimp][1]=CandiList[0][1];
				CandiList[0][0]=i;
				CandiList[0][1]=j;
				copywc(newwc,bestwc,num_wc);
				nimp++;
			}
			else {
				CandiList[nimp][0]=i;
				CandiList[nimp][1]=j;
				nimp++;
			}
		}
	}
}

void newfill_CL(int *oldwc, int *newwc, int *bestwc, int num_wc, int &numI, int &numN, int i, int j){
	//	double judge_f;
		int judge_f;
		judge_f=judge_part(oldwc,newwc,num_wc);
		if(judge_f<0){								//strictly better than the previous one
			if(numI==0){
				CandiList[0][0]=i;
				CandiList[0][1]=j;
				copywc(newwc,bestwc,num_wc);
				numI++;
			}
			else {
				if(judge_part(bestwc,newwc,num_wc)<0){
					CandiList[numI][0]=CandiList[0][0];
					CandiList[numI][1]=CandiList[0][1];
					CandiList[0][0]=i;
					CandiList[0][1]=j;
					copywc(newwc,bestwc,num_wc);
					numI++;
				}
				else {
					CandiList[numI][0]=i;
					CandiList[numI][1]=j;
					numI++;
				}
			}
		}
		else {
			Candidatemove[numN][0]=i;
			Candidatemove[numN][1]=j;
			numN++;
		}
}

void balance_CL(int i, int j){

	CandiList[0][0]=i;
	CandiList[0][1]=j;
}

void copywc(int *souwc, int *deswc, int l){
	int i=0;
	for(i=0;i<l;i++) deswc[i]=souwc[i];
}


int calcul_tenure(int ite){
	int aj;
	int tenure;
	int index_j;
	index_j=floor((ite%1500)/100);
	aj=list_aj[index_j];
	tenure=(aj-1)*100+rand()%100;
	return tenure;
}

void update_TL(int node1, int node2, int ite){
	TL[node1][node2]=ite+calcul_tenure(ite);
	TL[node2][node1]=ite+calcul_tenure(ite);
}

int RepairNF1(int &ite, Struc_Sol *C, Struc_Sol *LO, int th){
	int nimp=0;
	int i=0,j;
	int cb=C->cbmp;
	int label_midu=-1;
	int dis_u=-1;
	int set_s[num_vex];
	int taille_s=0;
	int max_wc=num_vex/2;
	int num_wc=num_vex/2+1;
	int newwc[num_wc];
	int oldwc[num_wc];
	int bestwc[num_wc];
	int u,v;
//	int oldcbs;
	int newcbs;
//	int bestcbs;

//	oldcbs=C->cbs;

	//decide the critical node
	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]>=cb) CV[i]=1;
		else CV[i]=0;

	}
	for(i=0;i<num_wc;i++) oldwc[i]=C->wc[i];

	if(cb+th<max_wc) oldwc[cb+th]=1;
	else oldwc[max_wc]=num_vex;

	for(i=0;i<num_vex;i++){
		if(CV[i]==0) continue;
		CV[i]=0;
		u=i;
		label_midu=find_midu(C,u);
		dis_u=get_CycD(label_midu,C->permutation[u]);
		//to get the neighbor for i
		taille_s=0;
		for(j=0;j<dis_u;j++){
			if((label_midu+j+num_vex)%num_vex!=C->permutation[u]) {
				set_s[taille_s]=C->permutationNew[(label_midu+j+num_vex)%num_vex];
				taille_s++;
			}
			if((label_midu-j+num_vex)%num_vex!=C->permutation[u]) {
				set_s[taille_s]=C->permutationNew[(label_midu-j+num_vex)%num_vex];
				taille_s++;
			}
		}
		//browse all the node for i
		for(j=0;j<taille_s;j++){
			v=set_s[j];
			if (TL[u][v]>=ite) continue;
			if (TL[v][u]>=ite) continue;

			get_newwc(C,newwc,u,v,newcbs);
//			cbs_Candidate(oldwc,newwc,bestwc,num_wc,oldcbs,newcbs,bestcbs,nimp,u,v);
			fill_CL(oldwc,newwc,bestwc,num_wc,nimp,u,v);
		}
	}
	return nimp;
}

int RepairNF2(int &ite, Struc_Sol *C, Struc_Sol *LO, int th){
	int nimp=0;
	int i=0,j;
	int cb=C->cbmp;

	int set_s[num_vex];

	int num_wc=num_vex/2+1;
	int newwc[num_wc];
	int oldwc[num_wc];
	int bestwc[num_wc];
	int u,v;
	int newcbs;


	int taille_s=num_vex*percent_can;
	//decide the critical node
	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]>=cb) CV[i]=1;
		else CV[i]=0;
	}
	for(i=0;i<num_wc;i++) oldwc[i]=C->wc[i];
	if(cb+th<num_vex/2) oldwc[cb+th]=1;
	else oldwc[num_vex]=1;

	for(i=0;i<num_vex;i++){
		if(CV[i]==0) continue;
		CV[i]=0;
		u=i;
		//generate the exc for each candidate
		generate_sets(set_s,num_vex);
		//browse all the node for i
		for(j=0;j<taille_s;j++){
			v=set_s[j];
			if(u==v) continue;
			if (TL[u][v]>=ite) continue;
			if (TL[v][u]>=ite) continue;

			get_newwc(C,newwc,u,v,newcbs);
//			cbs_Candidate(oldwc,newwc,bestwc,num_wc,oldcbs,newcbs,bestcbs,nimp,u,v);
			fill_CL(oldwc,newwc,bestwc,num_wc,nimp,u,v);
		}
	}
	return nimp;
}

int RepairNF3(int & ite, Struc_Sol *C, Struc_Sol *LO, int th){
	return 0;
}
int RepairNF4(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int &retbz){
	int nimp=0;
	int i=0,j;
	int cb=C->cbmp;

	int set_s[num_vex];

	int num_wc=num_vex/2+1;
	int newwc[num_wc];
//	int oldwc[num_wc];
	int bestwc[num_wc];
	int u,v;
	int bz=0;
	int bestbz=0;


	int taille_s=num_vex*percent_can;

	//decide the critical node
	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]>=cb*0.5) CV[i]=1;
		else CV[i]=0;
	}

	/*initial the wc*/
	for(i=0;i<num_wc;i++) bestwc[i]=num_edge;


	/*begin to find*/
	for(i=0;i<num_vex;i++){
		if(CV[i]==0) continue;
		CV[i]=0;
		u=i;
		//generate the exc for each candidate
		generate_sets(set_s,num_vex);
		//browse all the node for i
		for(j=0;j<taille_s;j++){
			v=set_s[j];
			if (TL[u][v]>=ite) continue;
			if (TL[v][u]>=ite) continue;

			bz=find_balance(C,newwc,u,v);
			if(bz>=0) continue;
			nimp=1;
			if(judge_part(bestwc,newwc,num_wc)<=0){
				balance_CL(u,v);
				bestbz=bz;
				for(int k=0;k<num_wc;k++) bestwc[k]=newwc[k];
			}
		}
	}
	retbz=bestbz;
	return nimp;
}


int RepairNF5(int & ite, Struc_Sol *C, Struc_Sol *LO, int th){
	int nimp=0;
	int i=0,j;
	int cb=C->cbmp;
	int label_midu=-1;
	int dis_u=-1;
	int set_s[num_vex];
	int taille_s=0;
	int num_wc=num_vex/2+1;
	int newwc[num_wc];
	int oldwc[num_wc];
	int bestwc[num_wc];
	int u,v;
	int oldcbs;
	int newcbs;
	int bestcbs;

	oldcbs=C->cbs;

	//decide the critical node
	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]>=cb) CV[i]=1;
		else CV[i]=0;
		oldwc[i]=C->wc[i];
	}
	if(cb+th<num_vex/2) oldwc[cb+th]=1;
	else oldwc[num_vex]=1;

	for(i=0;i<num_vex;i++){
		if(CV[i]==0) continue;
		CV[i]=0;
		u=i;
		label_midu=find_midu(C,u);
		dis_u=get_CycD(label_midu,C->permutation[u]);
		//to get the neighbor for i
		taille_s=0;
		for(j=0;j<=dis_u;j++){
			if((label_midu+j+num_vex)%num_vex!=C->permutation[u]) {
				set_s[taille_s]=C->permutationNew[(label_midu+j+num_vex)%num_vex];
				taille_s++;
			}
			if((label_midu-j+num_vex)%num_vex!=C->permutation[u]) {
				set_s[taille_s]=C->permutationNew[(label_midu-j+num_vex)%num_vex];
				taille_s++;
			}
		}
		//browse all the node for i
		for(j=0;j<taille_s;j++){
			v=set_s[j];
			if (u==v) continue;
			if (TL[u][v]>=ite) continue;
			if (TL[v][u]>=ite) continue;

			get_newwc(C,newwc,u,v,newcbs);
			cbs_Candidate(oldwc,newwc,bestwc,num_wc,oldcbs,newcbs,bestcbs,nimp,u,v);
//			fill_CL(oldwc,newwc,bestwc,num_wc,nimp,u,v);
		}
	}
	return nimp;
}
int RepairNF6(int & ite, Struc_Sol *C, Struc_Sol *LO, int th){
	int i,j,k;
	int cb=C->cbmp;
	int num_wc=num_vex/2+1;
//	int newwc[num_wc];
//	int oldwc[num_wc];
	int bestwc[num_wc];
	int destinodes[num_vex];
	int N_des;
	int u,v;
	int sss,eee;
	int nimp=0;

	/*decide the critical nodes as the startnodes*/
	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]>=cb) CV[i]=1;
		else CV[i]=0;
	}
	for(i=0;i<num_wc;i++) bestwc[i]=C->wc[i];
	/*******************************************/
	CandiSI[0][0]=-1;
	CandiSI[0][1]=-1;
	CandiSI[0][2]=-1;
	/************************************/
	for(i=0;i<num_vex;i++){
		/*decide the destinodes for each startnodes*/
		if(CV[i]) u=i;
		else continue;
		N_des=0;
		sss=start_end[u][0];
		eee=start_end[u][1];
		for(j=sss;j<=eee;j++){
			if(get_CycD(C->permutation[u],C->permutation[edge[j][1]])==cb){
				destinodes[N_des]=edge[j][1];
				N_des++;
			}
		}
		if(N_des==0){
			cout<<"error in RepairNF6, calcul CycD error"<<endl;
			exit(5);
		}
		/********************************************/
		/*copy the current solution to the tempCS and operate the Shift-insert save the best*/
		for(j=0;j<N_des;j++){
			CopySolution(C,CurrentS);
			v=destinodes[j];
			for(k=0;k<cb-1;k++){
				SI_move(CurrentS,u,v,1,ite);
				if(judge_part(bestwc, CurrentS->wc,num_wc)<0) {
					copywc(CurrentS->wc, bestwc, num_wc);
					nimp++;
					CandiSI[0][0]=u;
					CandiSI[0][1]=v;
					CandiSI[0][2]=k+1;
				}
			}

		}
	}

	/*operate the SI move on the current solution*/
	return nimp;
}

int RepairMN1(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int & numI, int & numN){
//	int nimp=0;
	int i=0,j;
	int cb=C->cbmp;
	int label_midu=-1;
	int dis_u=-1;
	int set_s[num_vex];
	int taille_s=0;
//	int max_wc=num_vex/2;
	int num_wc=num_vex/2+1;
	int newwc[num_wc];
	int oldwc[num_wc];
	int bestwc[num_wc];
	int u,v;
//	int oldcbs;
	int newcbs;
//	int bestcbs;

//	oldcbs=C->cbs;
	numI=0;
	numN=0;

	//decide the critical node
	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]>=cb) CV[i]=1;
		else CV[i]=0;

	}
	for(i=0;i<num_wc;i++) oldwc[i]=C->wc[i];

//	if(cb+th<max_wc) oldwc[cb+th]=1;
//	else oldwc[max_wc]=num_vex;

	for(i=0;i<num_vex;i++){
		if(CV[i]==0) continue;
		CV[i]=0;
		u=i;
		label_midu=find_midu(C,u);
		dis_u=get_CycD(label_midu,C->permutation[u]);
		//to get the neighbor for i
		taille_s=0;
		for(j=0;j<dis_u;j++){
			if((label_midu+j+num_vex)%num_vex!=C->permutation[u]) {
				set_s[taille_s]=C->permutationNew[(label_midu+j+num_vex)%num_vex];
				taille_s++;
			}
			if((label_midu-j+num_vex)%num_vex!=C->permutation[u]) {
				set_s[taille_s]=C->permutationNew[(label_midu-j+num_vex)%num_vex];
				taille_s++;
			}
		}
		//browse all the node for i
		for(j=0;j<taille_s;j++){
			v=set_s[j];
			if (TL[u][v]>=ite) continue;
			if (TL[v][u]>=ite) continue;

			get_newwc(C,newwc,u,v,newcbs);
//			cbs_Candidate(oldwc,newwc,bestwc,num_wc,oldcbs,newcbs,bestcbs,nimp,u,v);
//			fill_CL(oldwc,newwc,bestwc,num_wc,nimp,u,v);
			newfill_CL(oldwc,newwc,bestwc,num_wc,numI,numN,u,v);
		}
	}
	return 1;
}
int RepairMN2(int & ite, Struc_Sol *C, Struc_Sol *LO, int th, int & numI, int & numN){
//	int nimp=0;
	int i=0,j;
	int cb=C->cbmp;

	int set_s[num_vex];

	int num_wc=num_vex/2+1;
	int newwc[num_wc];
	int oldwc[num_wc];
	int bestwc[num_wc];
	int u,v;
	int newcbs;
//	int oldcbs;
//	int bestcbs;

	int taille_s=num_vex*percent_can;
//	oldcbs=C->cbs;
	numI=0;
	numN=0;
	//decide the critical node
	for(i=0;i<num_vex;i++){
		if(C->cbnodes[i]>=cb) CV[i]=1;
		else CV[i]=0;
	}
	for(i=0;i<num_wc;i++) oldwc[i]=C->wc[i];

//	if(cb+th<num_wc-1) oldwc[cb+th]=1;
//	else oldwc[num_wc-1]=num_edge;

	for(i=0;i<num_vex;i++){
		if(CV[i]==0) continue;
		CV[i]=0;
		u=i;
		//generate the exc for each candidate
		generate_sets(set_s,num_vex);
		//browse all the node for i
		for(j=0;j<taille_s;j++){
			v=set_s[j];
			if(u==v) continue;
			if (TL[u][v]>=ite) continue;
			if (TL[v][u]>=ite) continue;

			get_newwc(C,newwc,u,v,newcbs);
//			cbs_Candidate(oldwc,newwc,bestwc,num_wc,oldcbs,newcbs,bestcbs,nimp,u,v);
//			fill_CL(oldwc,newwc,bestwc,num_wc,nimp,u,v);
			newfill_CL(oldwc,newwc,bestwc,num_wc,numI,numN,u,v);
		}
	}
	return 2;
}

int shift_insert(Struc_Sol *C){
	return 0;
}

int get_bestCB_pop(int value_alb, int num_pop){
	int i;
	int bestreturn=num_vex/2+1;

	for(i=0;i<num_pop;i++) if(PopS[i].cbmp<bestreturn) bestreturn=PopS[i].cbmp;

	return (bestreturn);
}
double get_entropy_pop(int value_alb, int num_pop){

	int i,j,k;

	double entropy=0.0;
	for(i=0;i<num_vex;i++) for(j=0;j<num_vex;j++) mat_sum[i][j]=0;
	for(i=0;i<num_vex;i++)
		for(k=0;k<num_pop;k++)
			mat_sum[i][PopS[k].permutation[i]]++;
	for(i=0;i<num_vex;i++) for(j=0;j<num_vex;j++) {
		if(mat_sum[i][j]==0) continue;
		entropy=entropy+1.0*mat_sum[i][j]/num_pop*log(mat_sum[i][j]*1.0/num_pop);
	}
	entropy=-entropy/(num_vex*log(num_vex));

	return entropy;
}
double get_averagedistance_pop(int value_alb, int num_pop){
	int i,j;
	double sum_dis=0.0;
	double ave_dis;
	for(i=0;i<num_pop;i++) for(j=i+1;j<num_pop;j++){
		sum_dis=sum_dis+get_dis_tsp(PopS[i].permutation,PopS[j].permutation,num_vex, mat_sum);
	}
	ave_dis=2*sum_dis/((num_pop-1)*num_pop);
	return ave_dis;
}

void IniSol(Struc_Sol *Sol){
	int i,temp;
	int expermu[num_vex];

	for (i=0;i<num_vex;i++) Sol->permutation[i]=i;
	for (i=0;i<num_vex;i++) expermu[i]=rand()%num_vex;
	for (i=0;i<num_vex;i++){
		temp=Sol->permutation[i];
		Sol->permutation[i]=Sol->permutation[expermu[i]];
		Sol->permutation[expermu[i]]=temp;
	}
	for (i=0;i<num_vex;i++) Sol->permutationNew[Sol->permutation[i]]=i;
	EvaSol(Sol);
}

void EvaSol(Struc_Sol *Sol){
	int ub_cb=num_vex/2+1;
	int i;
	int cd;

	Sol->cbs=0;
	Sol->cbmp=0;
	memset(Sol->wc,0,ub_cb*sizeof(int));
	memset(Sol->cbnodes,0,num_vex*sizeof(int));

	for (i=0;i<num_vex;i++) Sol->permutationNew[Sol->permutation[i]]=i;
	for (i=0;i<num_edge;i++){
		cd=get_CycD(Sol->permutation[simple_edge[i][0]],Sol->permutation[simple_edge[i][1]]);
		Sol->wc[cd]++;
		Sol->cbs=Sol->cbs+cd;													// the sum of cbs
		if(cd>Sol->cbnodes[simple_edge[i][0]]) Sol->cbnodes[simple_edge[i][0]]=cd;
		if(cd>Sol->cbnodes[simple_edge[i][1]]) Sol->cbnodes[simple_edge[i][1]]=cd;
		if(cd>Sol->cbmp) Sol->cbmp=cd;
	}
}



void CopySolution(Struc_Sol *SouSol, Struc_Sol *DesSol){
	memcpy(DesSol->permutation, SouSol->permutation, num_vex*sizeof(int));
	memcpy(DesSol->permutationNew, SouSol->permutationNew, num_vex*sizeof(int));
	memcpy(DesSol->wc, SouSol->wc, (num_vex/2+1)*sizeof(int));
	memcpy(DesSol->cbnodes, SouSol->cbnodes, num_vex*sizeof(int));
	DesSol->cbmp = SouSol->cbmp;
	DesSol->cbs=SouSol->cbs;
}

int get_CycD(int label1, int label2){
	int absD;
	int cycD;
	absD=abs(label1-label2);
	cycD=min(absD, num_vex-absD);
	return cycD;
}

void get_newwc(Struc_Sol *C, int *wc, int node1, int node2, int &ncbs){
	int lab_oldnode1,lab_oldnode2;
	int lab_newnode1,lab_newnode2;
	int i;
	int ls,le;
	int cycd;

	ncbs=C->cbs;									//get the old cbs at the beginning
	for(i=0;i<num_vex/2+1;i++) wc[i]=C->wc[i];
	lab_newnode2=lab_oldnode1=C->permutation[node1];
	lab_newnode1=lab_oldnode2=C->permutation[node2];
	/**************************************************************/
	ls=start_end[node1][0];
	le=start_end[node1][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode1!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode1, C->permutation[edge[i][1]]);
			wc[cycd]--;
			ncbs=ncbs-cycd;							//udpate cbs
			/*get the new cycd*/
			cycd=get_CycD(lab_newnode1, C->permutation[edge[i][1]]);
			wc[cycd]++;
			ncbs=ncbs+cycd;
		}
	}
	/**************************************************************/
	ls=start_end[node2][0];
	le=start_end[node2][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode2!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode2, C->permutation[edge[i][1]]);
			wc[cycd]--;
			ncbs=ncbs-cycd;
			/*get the new cycd*/
			cycd=get_CycD(lab_newnode2, C->permutation[edge[i][1]]);
			wc[cycd]++;
			ncbs=ncbs+cycd;
		}
	}
	/**************************************************************/
}

int find_balance(Struc_Sol *C,int *wc,int node1, int node2){
	int lab_oldnode1,lab_oldnode2;
	int lab_newnode1,lab_newnode2;
	int i;
	int ls,le;
	int cycdold;
	int cycdnew;
	int bz=0;
	int fazhi=alreadybest;


	lab_newnode2=lab_oldnode1=C->permutation[node1];
	lab_newnode1=lab_oldnode2=C->permutation[node2];
	/**************************************************************/
	ls=start_end[node1][0];
	le=start_end[node1][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode1!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycdold=get_CycD(lab_oldnode1, C->permutation[edge[i][1]]);
			wc[cycdold]--;
			/*get the new cycd*/
			cycdnew=get_CycD(lab_newnode1, C->permutation[edge[i][1]]);
			wc[cycdnew]++;
			/*update the bz*/
			if((cycdold<=fazhi)&&(cycdnew>fazhi)) bz++;
			else if((cycdold>fazhi)&&(cycdnew<=fazhi)) bz--;

		}
	}
	/**************************************************************/
	ls=start_end[node2][0];
	le=start_end[node2][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode2!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycdold=get_CycD(lab_oldnode2, C->permutation[edge[i][1]]);
			wc[cycdold]--;
			/*get the new cycd*/
			cycdnew=get_CycD(lab_newnode2, C->permutation[edge[i][1]]);
			wc[cycdnew]++;
			/*update the bz*/
			if((cycdold<=fazhi)&&(cycdnew>fazhi)) bz++;
			else if((cycdold>fazhi)&&(cycdnew<=fazhi)) bz--;

		}
	}
	/**************************************************************/
	return bz;
}



void get_swap(Struc_Sol *C, int node1, int node2){
	int temp_node;
	int i;
	int newwc[num_vex/2+1];
	//update the newwc and cbnodes
	update_swap(C, node1, node2, newwc);

	//update the wc and cbmp
	for (i=0;i<num_vex/2+1;i++){
		C->wc[i]=newwc[i];
	}
	for (i=num_vex/2;newwc[i]==0;i--);
	C->cbmp=i;
	//update the permutation and permutationNew
	temp_node=C->permutation[node1];
	C->permutation[node1]=C->permutation[node2];
	C->permutation[node2]=temp_node;

	C->permutationNew[C->permutation[node1]]=node1;
	C->permutationNew[C->permutation[node2]]=node2;

	update_cbnodes(C,node1,node2);
}

void update_swap(Struc_Sol *C, int node1, int node2,int *wc){
	int lab_oldnode1,lab_oldnode2;
	int lab_newnode1,lab_newnode2;
	int i;
	int ls,le;
	int cycd;

	for(i=0;i<num_vex/2+1;i++) wc[i]=C->wc[i];
	lab_newnode2=lab_oldnode1=C->permutation[node1];
	lab_newnode1=lab_oldnode2=C->permutation[node2];
	/**************************************************************/
	ls=start_end[node1][0];
	le=start_end[node1][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode1!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode1, C->permutation[edge[i][1]]);
			wc[cycd]--;
			C->cbs=C->cbs-cycd;					//delete the old edge

			/*get the new cycd*/
			cycd=get_CycD(lab_newnode1, C->permutation[edge[i][1]]);
			wc[cycd]++;
			C->cbs=C->cbs+cycd;					//add the new edge
		}
	}
	/**************************************************************/
	ls=start_end[node2][0];
	le=start_end[node2][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode2!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode2, C->permutation[edge[i][1]]);
			wc[cycd]--;
			C->cbs=C->cbs-cycd;					//delete the old edge

			/*get the new cycd*/
			cycd=get_CycD(lab_newnode2, C->permutation[edge[i][1]]);
			wc[cycd]++;
			C->cbs=C->cbs+cycd;					//add the new edge
		}
	}
	/**************************************************************/
}

void update_cbnodes(Struc_Sol *C, int node1, int node2){
	int i,j,cbNi,cbNj,vl,ul;
//	int temp_node;

	cbNi=0;

	for(i=start_end[node1][0];i<=start_end[node1][1];i++){
		cbNj=0;
		for(j=start_end[edge[i][1]][0];j<=start_end[edge[i][1]][1];j++){
			vl=get_CycD(C->permutation[edge[i][1]], C->permutation[edge[j][1]]);
			if(vl>cbNj) cbNj=vl;
		}
		C->cbnodes[edge[i][1]]=cbNj;
		vl=get_CycD(C->permutation[node1], C->permutation[edge[i][1]]);
		if(vl>cbNi) cbNi=vl;
	}
	C->cbnodes[node1]=cbNi;

	cbNi=0;
	for(i=start_end[node2][0];i<=start_end[node2][1];i++){
		cbNj=0;
		for(j=start_end[edge[i][1]][0];j<=start_end[edge[i][1]][1];j++){
			ul=get_CycD(C->permutation[edge[i][1]], C->permutation[edge[j][1]]);
			if(ul>cbNj) cbNj=ul;
		}
		C->cbnodes[edge[i][1]]=cbNj;
		ul=get_CycD(C->permutation[node2], C->permutation[edge[i][1]]);
		if(ul>cbNi) cbNi=ul;
	}
	C->cbnodes[node2]=cbNi;

}

int find_midu(Struc_Sol *s, int u){

	int j,k,l;
//	int t;
	int taille_nei,flag_s,flag_e,label_insert,node_insert;
	int ne_nodes[num_vex];
//	int set_s[num_vex];
	int lu,ru,midu,dif,dif_max;
//	int nodec;
//	int taille_s;

	/*Create sequence of label of all the neighborhood verties of u*/
	taille_nei=0;
	flag_s=start_end[u][0];
	flag_e=start_end[u][1];
	for (j=flag_s;j<=flag_e;j++){
		node_insert=edge[j][1];
		label_insert=s->permutation[node_insert];
		if (taille_nei==0){
			ne_nodes[taille_nei]=label_insert;
		}
		else {
			for (k=0;k<taille_nei;k++){
				if (k==taille_nei-1){
					if (label_insert>ne_nodes[k]){
						ne_nodes[k+1]=label_insert;
						break;
					}
					else{
						ne_nodes[k+1]=ne_nodes[k];
						ne_nodes[k]=label_insert;
						break;
					}
				}
				else{
					if (label_insert>ne_nodes[k]){
						continue;
					}
					else{
						for (l=taille_nei-1;l>=k;l--){
							ne_nodes[l+1]=ne_nodes[l];
						}
							ne_nodes[k]=label_insert;
							break;
					}
				}
			}
		}
		taille_nei++;
	}
	ne_nodes[taille_nei]=ne_nodes[0]+num_vex;
	taille_nei++;
	/*look for the lu and ru, then calcul the midu and then looke for the vertex closer to midu*/
	lu=0;
	ru=0;
	midu=0;
//	nodec=0;
	dif=0;
	dif_max=0;
	for (j=0;j<taille_nei-1;j++){
		dif=ne_nodes[j+1]-ne_nodes[j];
		if (dif>dif_max){
			lu=ne_nodes[j+1];
			ru=ne_nodes[j];
			dif_max=dif;
		}
	}
	midu=(lu+ru+num_vex)/2;
	midu=midu%num_vex;
	return midu;
}


