//============================================================================
// Name        : BCP_mem.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
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
#include "BCP_X.h"
//#include "BCP_struct.h"
#include "BCP_c_LS.h"

using namespace std;
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

//extern function
extern int *get_vector(int size);
extern int **get_matrix(int num_row,int num_col);
extern int get_CycD(int leabel1, int label2);
//extern variable
extern int num_vex;
//extern int num_unk;
extern int num_edge;
//extern int max_deg;
extern int **start_end;
extern int **edge;
//extern int *degree_node;
//extern int **simple_edge;
extern char name_final_result[256];
//extern char file_G[256];
//extern char file_E[256];
//extern char file_D[256];
extern int alreadybest;

extern int L1;
extern int L2;
extern int L3;
//extern double alpha;

/***********************/

extern int maxi_cbs;

Struc_Sol *CurrentS,*BestS, *PopS, *ChS1, *ChS2;
Struc_Sol *LBS;
Struc_Sol *lb_local;
Struc_Sol *llb;
Struc_Sol *petlb;

// Variables time
struct tms glo_start;
struct tms glo_end;
struct tms midt;
double clockTicksPerSecond;
double startTimeSeconds;
double cpuTimeTheBest;

int **CandiSI;

//Local variable
int *CV;
int **Candidatemove;
int **mat_sum;
int **TL;
int **CandiList;
int **mat_dis;
int cuttingtime;

//declaration of function
int memSearch(int value_alb, int num_pop);
void setdataStructSol(int pop_sol);
void freedataStructSol(int pop_sol);

double get_time(void);

void generate_child(int num_pop,Struc_Sol *C1,Struc_Sol *C2, int length);
void quality_update_pop(Struc_Sol *C, int length);

int Shift_I(Struc_Sol *C, int &ite);
int Compromise_phase(Struc_Sol *C, Struc_Sol* LB, int &ite);
void ini_LS(Struc_Sol *souS, Struc_Sol *desS);
void SSearch(Struc_Sol *C);
int descent_upgrade(Struc_Sol *C, Struc_Sol* LB, int &ite, int th);
int descent_biggraph(Struc_Sol *C, Struc_Sol* LB, int &ite, int th);

void update_best(Struc_Sol *C, int ite);
int time_check(void);
void DistanceControlUpdatePopulation(Struc_Sol *C, int length,int pop_sol);
void ini_matdis(int pop_sol);
double calcul_Aslash(double ymax, double ymin, double y);
//local function

int RepairBig(int & ite, Struc_Sol *C,  int th);
Struc_Sol *GeneStrucSol(void);
void FreeStrucSol(Struc_Sol *poi);
int checkbest(void);

// Function
int checkbest(void){
	if(BestS->cbmp<=alreadybest) return 1;
	else return 0;
}

void FreeStrucSol(Struc_Sol *poi){
	free(poi->permutation);		poi->permutation=NULL;
	free(poi->permutationNew);	poi->permutationNew=NULL;
	free(poi->wc);				poi->wc=NULL;
	free(poi->cbnodes);			poi->cbnodes=NULL;
	free(poi);					poi=NULL;
}

Struc_Sol *GeneStrucSol(void){
	Struc_Sol *poi;
	int i;

	poi=(Struc_Sol*)malloc(sizeof(Struc_Sol));
	if(poi==NULL){
		cout<<"Memory error in poi"<<endl;
		exit(-1);
	}
	poi->permutation=(int*)get_vector(num_vex);
	poi->permutationNew=(int*)get_vector(num_vex);
	poi->wc=(int*)get_vector(num_vex/2+1);
	poi->cbnodes=(int*)get_vector(num_vex);

	for(i=0;i<num_vex;i++) {
		poi->permutation[i]=-1;
		poi->permutationNew[i]=-1;
		poi->cbnodes[i]=-1;
	}
	for(i=0;i<num_vex/2+1;i++) poi->wc[i]=-1;
	poi->cbmp=num_vex;
	poi->cbs=num_vex*num_vex;

	return poi;

}

int RepairBig(int & ite, Struc_Sol *C,  int th){
	int i,j;
	int num_wc=num_vex/2+1;
	int	precb;
	int prewc[num_wc];
	int u,v;
	int numI,numN;
	int temp;
	/***********************/
	int cb=C->cbmp;
	int label_midu=-1;
	int dis_u=-1;
	int set_s[num_vex];
	int taille_s=0;
	int newwc[num_wc];
	int oldwc[num_wc];
	int bestwc[num_wc];
//	int oldcbs;
	int newcbs;
//	int bestcbs;
	/**************************/

	CopySolution(C,llb);
	CopySolution(C,petlb);
	do{

		if(judge_part(llb->wc,petlb->wc,num_wc)<0) CopySolution(petlb,llb);
		precb=C->cbmp;
		copywc(C->wc,prewc,num_wc);
		for(i=0;i<num_vex;i++){
			if(C->cbnodes[i]>=cb) CV[i]=1;
			else CV[i]=0;
		}
//		for(i=0;i<num_wc;i++) oldwc[i]=C->wc[i];
		for(i=0;i<num_wc;i++) oldwc[i]=num_edge;

		for(i=0;i<num_vex;i++){
			if(CV[i]==0) continue;
			if(time_check()) return 0;
			CV[i]=0;
			u=i;
			/*find the midu*/
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
			numI=0;
			numN=0;
			/*browse each voisin*/
			for(j=0;j<taille_s;j++){
				v=set_s[j];
				if (TL[u][v]>=ite) continue;
				if (TL[v][u]>=ite) continue;
				get_newwc(C,newwc,u,v,newcbs);
				newfill_CL(oldwc,newwc,bestwc,num_wc,numI,numN,u,v);
			}
			/*get move*/
			ite++;
			if(numI!=0){
				temp=rand()%numI;
				temp=0;
				get_swap(C,CandiList[temp][0],CandiList[temp][1]);
				update_TL(CandiList[temp][0],CandiList[temp][1],ite);
			}
			else if(numI==0 && numN!=0){
				temp=rand()%numN;
				get_swap(C,Candidatemove[temp][0],Candidatemove[temp][1]);
				update_TL(Candidatemove[temp][0],Candidatemove[temp][1],ite);
			}
			/*update the CV*/
			for(int k=0;k<num_vex;k++)	if(C->cbnodes[i]<precb && CV[i]==1) CV[i]=0;

			/*update the petlb*/
			if(judge_part(petlb->wc,C->wc,num_wc)<0) CopySolution(C,petlb);

			/*update the BestS and time check*/
			update_best(C,ite);
//			if(time_check()) return 0;
		}
	}while(judge_part(llb->wc,petlb->wc,num_wc)<0);
	return 0;
}


double calcul_Aslash(double ymax, double ymin, double y){
	return 1.0*(y-ymin)/(ymax-ymin+1);
}

void ini_matdis(int pop_sol){
	int i,j;
	for(i=0;i<pop_sol;i++)
		for(j=0;j<pop_sol;j++) {
			mat_dis[i][j]=get_dis_tsp(PopS[i].permutationNew,PopS[j].permutationNew,num_vex,mat_sum);
		}
}

void DistanceControlUpdatePopulation(Struc_Sol *C, int length,int pop_sol){
	int i;
	int j;
	int ob_max;
	int ob_min;
	int	dismax;
	int dismin;
	double temp_reci;
	double obreciprocal_max;
	double obreciprocal_min;
	double score_pop[pop_sol];
	int distoP[pop_sol];
	int index_minscore=0;
	double minscore;
	/*item of new*/
	double score_new;
	int disnew2P;

	/*get obmax and min*/
	ob_max=C->cbmp;
	ob_min=C->cbmp;
	for(i=0;i<pop_sol;i++) {
		if(PopS[i].cbmp>ob_max) ob_max=PopS[i].cbmp;
		if(PopS[i].cbmp<ob_min) ob_min=PopS[i].cbmp;
	}
	obreciprocal_max=1.0/ob_min;
	obreciprocal_min=1.0/ob_max;

	/*remplir the mat_dis*/
	for(i=0;i<pop_sol;i++) {
		mat_dis[pop_sol][i]=get_dis_tsp(PopS[i].permutationNew,C->permutationNew,num_vex,mat_sum);
		mat_dis[i][pop_sol]=mat_dis[pop_sol][i];
	}


	/*get distoP and dismax min*/
	for(i=0;i<pop_sol;i++) distoP[i]=num_vex+1;					//adding code
	for(i=0;i<pop_sol;i++){
		for(j=0;j<pop_sol+1;j++){
			if(i==j) continue;
			if(distoP[i]>mat_dis[i][j]) distoP[i]=mat_dis[i][j];
		}
	}
	disnew2P=mat_dis[pop_sol][0];
	for(i=0;i<pop_sol;i++) {
		if(disnew2P>mat_dis[pop_sol][i]) disnew2P=mat_dis[pop_sol][i];
	}

	dismax=disnew2P;
	dismin=disnew2P;
	for(i=0;i<pop_sol;i++) {
		if(dismax<distoP[i]) dismax=distoP[i];
		if(dismin>distoP[i]) dismin=distoP[i];
	}



	/*get the score*/
	for(i=0;i<pop_sol;i++){
		temp_reci=1.0/PopS[i].cbmp;
		score_pop[i]=0.6*calcul_Aslash(obreciprocal_max,obreciprocal_min,temp_reci)+0.4*calcul_Aslash(dismax,dismin,distoP[i]);
	}
	temp_reci=1.0/C->cbmp;
	score_new=0.6*calcul_Aslash(obreciprocal_max,obreciprocal_min,temp_reci)+0.4*calcul_Aslash(dismax,dismin,disnew2P);


	/*find the min score*/
	minscore=score_pop[0];
	index_minscore=0;
	for(i=0;i<pop_sol;i++){
		if(minscore>score_pop[i]) {
			minscore=score_pop[i];
			index_minscore=i;
		}
	}


	/*decide to update and remplacer the worst one*/
	double pro_ran=rand()/(RAND_MAX+1.0);
	Struc_Sol *temp;

	if(minscore<score_new || pro_ran<0.3){
		temp=&PopS[index_minscore];
		CopySolution(C,temp);
		for(i=0;i<pop_sol;i++) mat_dis[index_minscore][i]=mat_dis[pop_sol][i];
		mat_dis[index_minscore][index_minscore]=0;									//adding code

	}
}


int time_check(void){
	double time_total=0.0;
	time_total=get_time();
	if(time_total>cuttingtime) return 1;
	else return 0;
}

void update_best(Struc_Sol *C, int ite){
	double time_total=0.0;
	time_total=get_time();
	if(C->cbmp<BestS->cbmp){
		CopySolution(C,BestS);
		cout<<"ite="<<ite<<" "<<BestS->cbmp<<endl;
		ofstream caout(name_final_result,ios::out|ios::app);
		if (caout.is_open()){
			caout<<ite<<" ";
			caout<<BestS->cbmp<<" ";
			caout<<time_total<<" ";
			caout<<endl;
			caout.close();
		}
	}
}


int descent_upgrade(Struc_Sol *C, Struc_Sol* LB, int &ite, int th){
	int count=0;
	double ra_best=0.0;
	int num_candi __attribute__ ((unused));
	int temp=-1;
	double judge;
	int numI;
	int numN;
	int flag=0;
	CopySolution(C,lb_local);

	while(count<L1){

		if(time_check()) return 0;
		numI=0;
		numN=0;

//		if(flag==1)	num_candi=RepairMN1(ite,C,LB,0,numI,numN);
//		else Compromise_move(C,LB,ite);

		num_candi=RepairMN1(ite,C,LB,0,numI,numN);
		// choose the move
		ite++;
		if (numI!=0 && numN!=0) {											// numI!=0, numN!=0
			ra_best=rand()/(RAND_MAX+1.0);
			if(ra_best>0) {
				temp=rand()%numI;
				temp=0;
				get_swap(C,CandiList[temp][0],CandiList[temp][1]);
				update_TL(CandiList[temp][0],CandiList[temp][1],ite);
			}
			else {
				temp=rand()%numN;
				get_swap(C,Candidatemove[temp][0],Candidatemove[temp][1]);
				update_TL(Candidatemove[temp][0],Candidatemove[temp][1],ite);
			}
		}
		else if(numI!=0 && numN==0){										//numI!=0, numN==0 take only from CandiList
			temp=rand()%numI;
			temp=0;
			get_swap(C,CandiList[temp][0],CandiList[temp][1]);
			update_TL(CandiList[temp][0],CandiList[temp][1],ite);
		}

		else if(numI==0 &&numN!=0){											//numI=0, numN!=0 take only from candidatemove
			temp=rand()%numN;
			get_swap(C,Candidatemove[temp][0],Candidatemove[temp][1]);
			update_TL(Candidatemove[temp][0],Candidatemove[temp][1],ite);
		}
																			//numI==0, numN==0
		/*******choose the next Nx************/
		if(flag==1){
			if(numI!=0) flag=1;
			else flag=2;
		}
		else flag=1;

		/*update the count*/
		if(judge_part(lb_local->wc,C->wc,num_vex/2+1)<0){
			count=0;
			CopySolution(C,lb_local);
		}
		else count++;
		/*update the best and LB*/
		judge=judge_part(BestS->wc,C->wc,num_vex/2+1);
		update_best(C,ite);

		if(judge<0) CopySolution(C,BestS);
		if(judge_part(LB->wc,C->wc,num_vex/2+1)<0) CopySolution(C,LB);
		if(checkbest()) return 0;
	}
	return 0;
}

int descent_biggraph(Struc_Sol *C, Struc_Sol* LB, int &ite, int th){
	int count=0;
//	double ra_best=0.0;
//	int num_candi=0;
//	int temp=-1;
	int flag=0;
	int judge;
	int numcandi __attribute__ ((unused));
	CopySolution(C,lb_local);

	while(count<L1){

		if(time_check()) return 0;

		numcandi=RepairBig(ite,C,0);
		// if(flag==1)	numcandi=RepairBig(ite,C,0);
		// else Compromise_move(C,LB,ite);
																			//numI==0, numN==0
		/*******choose the next Nx************/
		if(flag==1) flag=2;
		else flag=1;

		/*update the count*/
		if(judge_part(lb_local->wc,llb->wc,num_vex/2+1)<0){					//llb is the local optimal in RepairBig
			count=0;
			CopySolution(llb,lb_local);
		}
		else count++;
		/*update the best and LB*/
		judge=judge_part(BestS->wc,C->wc,num_vex/2+1);
		update_best(C,ite);


		if(judge<0) CopySolution(C,BestS);
		if(judge_part(LB->wc,lb_local->wc,num_vex/2+1)<0) CopySolution(lb_local,LB);
		if(checkbest()) return 0;
	}
	return 0;
}

int Shift_I(Struc_Sol *C, int &ite){
	int i;
	int nn_cv=0;
	int tC[num_vex];
	int index_chosen;
	int node_chosen;
	int temp_node;
	int pair_node=-1;
	int ss,ee;

	for(i=0;i<num_vex;i++)
		if(C->cbnodes[i]>=C->cbmp){
			tC[nn_cv]=i;
			nn_cv++;
		}
	index_chosen=rand()%nn_cv;
	node_chosen=tC[index_chosen];
	ss=start_end[node_chosen][0];
	ee=start_end[node_chosen][1];

	for(i=ss;i<=ee;i++){
		temp_node=edge[i][1];
		if(get_CycD(C->permutation[node_chosen],C->permutation[temp_node])==C->cbmp) {
			pair_node=temp_node;
			break;
		}
	}

	if(pair_node==-1) {
		cout<<"error in pairnode"<<endl;
		exit(-1);
	}
	if(rand()/(RAND_MAX+1.0)>=0.5) SI_move(C,node_chosen,pair_node,C->cbmp-1,ite);
	else SI_move(C,pair_node,node_chosen,C->cbmp-1,ite);

	return 0;
}

int Compromise_phase(Struc_Sol *C, Struc_Sol* LB, int &ite){
	int count=0;

//	int num_candi=0;
//	int temp=-1;
//	int lmax=len_th*num_vex;
	double judge;
//	int temp_bz;

	while(count<L2){

		if(time_check()) return 0;

		Compromise_move(C,LB,ite);

		judge=judge_part(BestS->wc,C->wc,num_vex/2+1);
		update_best(C,ite);

		if(judge<0) CopySolution(C,BestS);
		if(judge_part(LB->wc,C->wc,num_vex/2+1)<0) CopySolution(C,LB);
		count++;
	}

	return 0;
}

void SSearch(Struc_Sol *C){
	int ite=1;
//	double time_total=0;
	int count=0;
//	int th=0;
//	int cc=0;

	ini_LS(C,LBS);
	if(time_check()) return;

//	while(cc<2){
		count=0;
		while(count<L3){

//			Descent(C,LBS,ite,0);
//			newDescent(C,LBS,ite,0);
//			modified_descent(C,LBS,ite,0);
//			descent_discuss(C,LBS,ite,0);
//			descent_upgrade(C,LBS,ite,0);
			descent_biggraph(C,LBS,ite,0);

			if(checkbest()) return;
			if(time_check()) return;
//			cout<<"after Decent="<<C->cbmp<<endl;
			Compromise_phase(C,LBS,ite);

//			th=ceil(BestS->cbmp*(1/(0.00891104*BestS->cbmp+0.52663736)+0.16331589));
//			Threshold(C,LBS,ite,th);
//			Balance(C,LBS,ite,0);
//			cout<<"after balance="<<C->cbmp<<endl;
			count++;
			if(time_check()) return;
		}

//		cc++;
//	}

}

void ini_LS(Struc_Sol *souS, Struc_Sol *desS){
	int i,j;
	for(i=0;i<num_vex;i++) for(j=0;j<num_vex;j++) TL[i][j]=0;
	CopySolution(souS,desS);
}

void setdataStructSol(int pop_sol){
	//initial *CV,*CurrentS, *BestS, **Candidatemove,*PopS, *ChS1, *ChS2;
//	int ub_cb=num_vex/2+1 __attribute__((unused));

	int i;
	CV=(int*)get_vector(num_vex);
	for(i=0;i<num_vex;i++) CV[i]=-1;
	/*************************************/
	CurrentS=(Struc_Sol*)GeneStrucSol();

	/*************************************/
	BestS=(Struc_Sol*)GeneStrucSol();

	/*************************************/
	Candidatemove=(int**)get_matrix(num_vex*num_vex,2);
	/*************************************/

	/*************************************/
	ChS1=(Struc_Sol*)GeneStrucSol();


	/*************************************/
	ChS2=(Struc_Sol*)GeneStrucSol();


	/**intitial the mat_s********************/
	mat_sum=(int**)get_matrix(num_vex,num_vex);
	/**intitial the TL and CandiList********************/
	TL=(int **)get_matrix(num_vex,num_vex);
	CandiList=(int **)get_matrix(num_vex*num_vex,2);

	/* initial LBS*/
	LBS=(Struc_Sol*)GeneStrucSol();

	/**************initial matrix_distance***********************/
	mat_dis=(int**)get_matrix(pop_sol+1,pop_sol+1);
	/***********intial the candilist of shift insert*******************/
	CandiSI=(int **)get_matrix(1,3);
	/**initial lb_local**/
	lb_local=(Struc_Sol*)GeneStrucSol();


	/**initial llb**/
	llb=(Struc_Sol*)GeneStrucSol();


	/**initial petlb**/
	petlb=(Struc_Sol*)GeneStrucSol();

}

void freedataStructSol(int pop_sol){
	int i;
	//free the space of *CV, *CurrentS, *BestS, **Candidatemove
	free(CV);		CV=NULL;

	free(CurrentS->permutation);		CurrentS->permutation=NULL;
	free(CurrentS->permutationNew);		CurrentS->permutationNew=NULL;
	free(CurrentS->wc);					CurrentS->wc=NULL;
	free(CurrentS->cbnodes);			CurrentS->cbnodes=NULL;
	free(CurrentS);						CurrentS=NULL;

	free(BestS->permutation);			BestS->permutation=NULL;
	free(BestS->permutationNew);		BestS->permutationNew=NULL;
	free(BestS->wc);					BestS->wc=NULL;
	free(BestS->cbnodes);				BestS->cbnodes=NULL;
	free(BestS);						BestS=NULL;

	for(i=0;i<num_vex*num_vex;i++){
		free(Candidatemove[i]);
		Candidatemove[i]=NULL;
	}
	free(Candidatemove);				Candidatemove=NULL;

//	for(i=0;i<pop_sol;i++){
//		free(PopS[i].permutation);			PopS[i].permutation=NULL;
//		free(PopS[i].permutationNew);		PopS[i].permutationNew=NULL;
//		free(PopS[i].wc);					PopS[i].wc=NULL;
//		free(PopS[i].cbnodes);				PopS[i].cbnodes=NULL;
//	}
//	free(PopS);							PopS=NULL;

	free(ChS1->permutation);			ChS1->permutation=NULL;
	free(ChS1->permutationNew);			ChS1->permutationNew=NULL;
	free(ChS1->wc);						ChS1->wc=NULL;
	free(ChS1->cbnodes);				ChS1->cbnodes=NULL;
	free(ChS1);							ChS1=NULL;

	free(ChS2->permutation);			ChS2->permutation=NULL;
	free(ChS2->permutationNew);			ChS2->permutationNew=NULL;
	free(ChS2->wc);						ChS2->wc=NULL;
	free(ChS2->cbnodes);				ChS2->cbnodes=NULL;
	free(ChS2);							ChS2=NULL;

	for(i=0;i<num_vex;i++){
		free(mat_sum[i]);
		mat_sum[i]=NULL;
	}
	free(mat_sum);						mat_sum=NULL;
	/*free the TL*/
	for(i=0;i<num_vex;i++){
		free(TL[i]);
		TL[i]=NULL;
	}
	free(TL);							TL=NULL;
	/*free the CandiList*/
	for(i=0;i<num_vex*num_vex;i++){
		free(CandiList[i]);
		CandiList[i]=NULL;
	}
	free(CandiList);					CandiList=NULL;
	/*free the LBS*/
	free(LBS->permutation);				LBS->permutation=NULL;
	free(LBS->permutationNew);			LBS->permutationNew=NULL;
	free(LBS->wc);						LBS->wc=NULL;
	free(LBS->cbnodes);					LBS->cbnodes=NULL;
	free(LBS);							LBS=NULL;
	/*free the mat_dis*/
	for(i=0;i<pop_sol+1;i++){
		free(mat_dis[i]);
		mat_dis[i]=NULL;
	}
	free(mat_dis);						mat_dis=NULL;
	/*free the CandiSI*/
	for(i=0;i<1;i++){
		free(CandiSI[i]);
		CandiSI[i]=NULL;
	}
	free(CandiSI);						CandiSI=NULL;
	/*free the lb_local*/
	free(lb_local->permutation);		lb_local->permutation=NULL;
	free(lb_local->permutationNew);		lb_local->permutationNew=NULL;
	free(lb_local->wc);					lb_local->wc=NULL;
	free(lb_local->cbnodes);			lb_local->cbnodes=NULL;
	free(lb_local);						lb_local=NULL;
	/*free the llb*/
	free(llb->permutation);				llb->permutation=NULL;
	free(llb->permutationNew);			llb->permutationNew=NULL;
	free(llb->wc);						llb->wc=NULL;
	free(llb->cbnodes);					llb->cbnodes=NULL;
	free(llb);							llb=NULL;
	/*free the petlb*/
	free(petlb->permutation);				petlb->permutation=NULL;
	free(petlb->permutationNew);			petlb->permutationNew=NULL;
	free(petlb->wc);						petlb->wc=NULL;
	free(petlb->cbnodes);					petlb->cbnodes=NULL;
	free(petlb);							petlb=NULL;
}

double get_time(void){
	double endTimeSeconds=0.0;
	times(&glo_end);
	endTimeSeconds = glo_end.tms_utime/clockTicksPerSecond;
	return (endTimeSeconds - startTimeSeconds);
}

void generate_child(int num_pop, Struc_Sol *C1, Struc_Sol *C2, int length){
	int idx_f,idx_m;
	Struc_Sol *tempSol __attribute((unused));
	Struc_Sol *tSolF, *tSolM;

	do{
		idx_f=rand()%num_pop;
		idx_m=rand()%num_pop;
	}while(idx_f==idx_m);
//	PMX(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew);
//	CX(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew);
//	OX(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew);
//	OX2(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew);
//	while(DPX(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew, mat_sum)==0){}
//	OX_modified(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew, PopS[idx_f].cbnodes, PopS[idx_m].cbnodes, PopS[idx_f].permutation, PopS[idx_m].permutation,PopS[idx_f].cbmp);
	tSolF=&PopS[idx_f];
	tSolM=&PopS[idx_m];
	CLONEReconstruct(tSolF,tSolM,C1);

	for (int i=0;i<num_vex;i++) C1->permutation[C1->permutationNew[i]]=i;
//	for (int i=0;i<num_vex;i++) C2->permutation[C2->permutationNew[i]]=i;

	EvaSol(C1);
//	cout<<"finish the generate children"<<endl;
//	EvaSol(C2);
}

void quality_update_pop(Struc_Sol *C,int length){
	int i;
	int index_pire=-1, max_dif=0;
	Struc_Sol *temp;

	for(i=0;i<length;i++) if(PopS[i].cbmp-C->cbmp>max_dif) {
		index_pire=i;
		max_dif=PopS[i].cbmp-C->cbmp;
	}
	if(index_pire!=-1){
		temp=&PopS[index_pire];
		CopySolution(C,temp);
	}
}

int memSearch(int value_alb, int num_pop){
	int i __attribute__ ((unused));
	int bc_final;
	double time_total=0;
	Struc_Sol *temp __attribute__ ((unused));
	int count=0;


	maxi_cbs=num_edge*num_vex/2;						// new maxi_cbs

    clockTicksPerSecond = (double)sysconf(_SC_CLK_TCK); // Get the clock ticks per second
	setdataStructSol(num_pop);

	cuttingtime=600;

	times(&glo_start);
	startTimeSeconds = glo_start.tms_utime/clockTicksPerSecond;

	IniSol(BestS);

	time_total=get_time();
	ofstream detout("detail.txt",ios::out|ios::trunc);

	ofstream caout(name_final_result,ios::out|ios::app);
	if (caout.is_open()){
		caout<<"0 ";
		caout<<BestS->cbmp<<" ";
		caout<<time_total<<" ";
		caout<<endl;
		caout.close();
	}
	
	// Initialise the population
//	for(i=0;i<num_pop;i++){
//		if(time_check()) break;
//
//		temp=&PopS[i];
//		IniSol(temp);
//		SSearch(temp);
//		if(BestS->cbmp<=value_alb) {
//			EvaSol(BestS);
//			bc_final=BestS->cbmp;
//			freedataStructSol(num_pop);
//			return bc_final;
//		}
//		cout<<"initial "<<i<<" item"<<endl;
//	}
//	ini_matdis(num_pop);

	IniSol(ChS1);
	cout<<"finish all initial"<<endl;
	//LS and memetic
	while(BestS->cbmp>value_alb){
		if(time_check()) break;
//		generate_child(num_pop, ChS1, ChS2, num_vex);
//		cout<<"crossover operated, the gene="<<count<<endl;
		SSearch(ChS1);
//		Modified_DFS(ChS1,ChS1);
		M2_DFS(ChS1,ChS1);
		EvaSol(ChS1);
//		DFSReconstruct(ChS1, ChS1);
//		EvaSol(ChS1);

//		quality_update_pop(LBS,num_pop);
//		DistanceControlUpdatePopulation(LBS, num_vex, num_pop);

		count++;
//		time_total=get_time();
//		if(time_total>=cuttingtime) break;
	}
	//verify the result and return the final cbmp
	EvaSol(BestS);
	bc_final=BestS->cbmp;
	freedataStructSol(num_pop);

	time_total=get_time();
	cout<<"endtime="<<time_total<<endl;
	return bc_final;
}
