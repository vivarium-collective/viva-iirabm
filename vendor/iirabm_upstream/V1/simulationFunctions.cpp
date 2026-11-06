#include <stdlib.h>
#include <vector>
#include <iostream>
#include "agents.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include "Parameters.h"
using namespace std;

extern mt19937 generator;
extern uniform_int_distribution<int> distribution100,distribution3;
extern int cellGrid[101][101],procID;
extern float system_oxy,oxyDeficit,totalInfection,total_TNF,total_sTNFr,total_IL10,
      	total_IL6,total_GCSF,total_proTH1,total_proTH2,oxyheal,total_IFNg,total_PAF,
      	total_IL1,total_IL4,total_IL8,total_IL12,total_sIL1r,total_IL1ra;
extern vector<EC> ecArray;
extern vector<pmn> pmnArray;
extern vector<mono> monoArray;
extern vector<TH0> TH0array;
extern vector<TH1> TH1array;
extern vector<TH2> TH2array;
extern vector<pmn_marrow> pmn_marrowArray;
extern vector<mono_marrow> mono_marrowArray;
extern vector<TH0_germ> TH0_germArray;
extern vector<TH1_germ> TH1_germArray;
extern vector<TH2_germ> TH2_germArray;

extern float PAFmult,TNFmult,sTNFrmult,IL1mult,sIL1rmult,IL1ramult,IFNgmult,IL4mult,
	IL8mult,IL10mult,IL12mult,GCSFmult;
extern float myIntervention[genomeLength],testInt[testLength];
extern float allInterventions[ga_popSize][genomeLength];

extern float mGCSF,mPAF,mTNF,mSTNFR,mIL1,mSIL1R,mIL1RA,mIFNg,mIL4,mIL8,mIL10,mIL12;

void wiggle(int* orientation, int* x, int* y);
void getAhead(int orient, int x, int y, int *xl, int *xm, int *xr, int *yl, int *ym, int *yr);
void adjustOrientation(int* orientation, int leftOrRight);
void move(int orient);
void updateSystemOxy(int step);
void evaporate();
void initialize();
void injure_infectionFRD(int inj_number);
void readIntervention();

void readParams(int procID, int nProc, int p1[], int* np1, float p2[], int* np2, int p3[], 
	int* np3, int p4[], int* np4,int p5[], int* np5, int p6[], int* np6);	

void diffuse(){
	float nFactor;
	int i,j,totalCells;
	
	totalCells=xDim*yDim;	
	nFactor=1.0/8.0;
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].endotoxin*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].endotoxin-=ecArray[i].endotoxin;
		ecArray[i].endotoxin+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].PAF*0.6*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].PAF-=(0.6*ecArray[i].PAF);
		ecArray[i].PAF+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].cytotox*0.4*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].cytotox-=(0.4*ecArray[i].cytotox);
		ecArray[i].cytotox+=ecArray[i].tempSignal;
	}

	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].TNF*0.6*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].TNF-=(0.6*ecArray[i].TNF);
		ecArray[i].TNF+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].sTNFr*0.8*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].sTNFr-=(0.8*ecArray[i].sTNFr);
		ecArray[i].sTNFr+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].IL1*0.6*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].IL1-=(0.6*ecArray[i].IL1);
		ecArray[i].IL1+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].IFNg*0.8*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].IFNg-=(0.8*ecArray[i].IFNg);
		ecArray[i].IFNg+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].IL8*0.6*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].IL8-=(0.6*ecArray[i].IL8);
		ecArray[i].IL8+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].IL10*0.8*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].IL10-=(0.8*ecArray[i].IL10);
		ecArray[i].IL10+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].IL1ra*0.8*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].IL1ra-=(0.8*ecArray[i].IL1ra);
		ecArray[i].IL1ra+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].sIL1r*0.8*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].sIL1r-=(0.8*ecArray[i].sIL1r);
		ecArray[i].sIL1r+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].IL12*0.8*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].IL12-=(0.8*ecArray[i].IL12);
		ecArray[i].IL12+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].IL4*0.8*nFactor);
		}
	}
	for(i=0;i<totalCells;i++){
		ecArray[i].IL4-=(0.8*ecArray[i].IL4);
		ecArray[i].IL4+=ecArray[i].tempSignal;
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].tempSignal=0;
//		cout<<"GCSF="<<ecArray[i].GCSF<<"\n";
	}
	for(i=0;i<totalCells;i++){
		for(j=0;j<8;j++){
			ecArray[ecArray[i].neighbors[j]].tempSignal+=(ecArray[i].GCSF*nFactor);	
		}
	}
	
	for(i=0;i<totalCells;i++){
		ecArray[i].GCSF-=ecArray[i].GCSF;
		ecArray[i].GCSF+=ecArray[i].tempSignal;
	}
}

void evaporate(){
	int size,i;
	size=ecArray.size();
	for(i=0;i<size;i++){
		ecArray[i].endotoxin*=0.7;
		if(ecArray[i].endotoxin<0.01){ecArray[i].endotoxin=0;}
		ecArray[i].PAF*=0.7;
		if(ecArray[i].PAF<0.01){ecArray[i].PAF=0;}
		ecArray[i].cytotox*=0.7;
		if(ecArray[i].cytotox<0.01){ecArray[i].cytotox=0;}
		ecArray[i].TNF*=0.8;
		if(ecArray[i].TNF<0.01){ecArray[i].TNF=0;}
		ecArray[i].IL1*=0.8;
		if(ecArray[i].IL1<0.01){ecArray[i].IL1=0;}
		ecArray[i].sTNFr*=0.9;
		if(ecArray[i].sTNFr<0.01){ecArray[i].sTNFr=0;}
		ecArray[i].IL1ra*=0.9;
		if(ecArray[i].IL1ra<0.01){ecArray[i].IL1ra=0;}
		ecArray[i].sIL1r*=0.9;
		if(ecArray[i].sIL1r<0.01){ecArray[i].sIL1r=0;}
		ecArray[i].IFNg*=0.8;
		if(ecArray[i].IFNg<0.01){ecArray[i].IFNg=0;}
		ecArray[i].IL8*=0.7;
		if(ecArray[i].IL8<0.01){ecArray[i].IL8=0;}
		ecArray[i].IL10*=0.95;
		if(ecArray[i].IL10<0.01){ecArray[i].IL10=0;}
		ecArray[i].IL12*=0.8;
		if(ecArray[i].IL12<0.01){ecArray[i].IL12=0;}
		ecArray[i].IL4*=0.95;
		if(ecArray[i].IL4<0.01){ecArray[i].IL4=0;}
		ecArray[i].GCSF*=0.95;
		if(ecArray[i].GCSF<0.01){ecArray[i].GCSF=0;}
	}
}

void updateSystemOxy(int step){
	int size,i;

	system_oxy=0;
	oxyDeficit=0;
	totalInfection=0;
	total_TNF=0;
	total_sTNFr=0;
	total_IL10=0;
	total_GCSF=0;
	total_proTH1=0;
	total_proTH2=0;
	total_IFNg=0;
	total_PAF=0;
	total_IL1=0;
	total_IL4=0;
	total_IL8=0;
	total_IL12=0;
	total_sIL1r=0;
	total_IL1ra=0;
	size=ecArray.size();
	for(i=0;i<size;i++){
		system_oxy+=(ecArray[i].oxy/100);
		totalInfection+=(ecArray[i].infection/100);
		total_TNF+=(ecArray[i].TNF/100);
		if(ecArray[i].TNF>mTNF){mTNF=ecArray[i].TNF;}
		total_sTNFr+=(ecArray[i].sTNFr/100);
		if(ecArray[i].sTNFr>mSTNFR){mSTNFR=ecArray[i].sTNFr;}
		total_IL10+=(ecArray[i].IL10/100);
		if(ecArray[i].IL10>mIL10){mIL10=ecArray[i].IL10;}
		total_GCSF+=(ecArray[i].GCSF/100);
		if(ecArray[i].GCSF>mGCSF){mGCSF=ecArray[i].GCSF;}
		total_IFNg+=(ecArray[i].IFNg/100);
		if(ecArray[i].IFNg>mIFNg){mIFNg=ecArray[i].IFNg;}
		total_PAF+=(ecArray[i].PAF/100);
		if(ecArray[i].PAF>mPAF){mPAF=ecArray[i].PAF;}
		total_IL1+=(ecArray[i].IL1/100);
		if(ecArray[i].IL1>mIL1){mIL1=ecArray[i].IL1;}
		total_IL4+=(ecArray[i].IL4/100);
		if(ecArray[i].IL4>mIL4){mIL4=ecArray[i].IL4;}
		total_IL8+=(ecArray[i].IL8/100);
		if(ecArray[i].IL8>mIL8){mIL8=ecArray[i].IL8;}
		total_IL12+=(ecArray[i].IL12/100);
		if(ecArray[i].IL12>mIL12){mIL12=ecArray[i].IL12;}
		total_sIL1r+=(ecArray[i].sIL1r/100);
		if(ecArray[i].sIL1r>mSIL1R){mSIL1R=ecArray[i].sIL1r;}
		total_IL1ra+=(ecArray[i].IL1ra/100);
		if(ecArray[i].IL1ra>mIL1RA){mIL1RA=ecArray[i].IL1ra;}
				
	}
	size=TH0array.size();
	for(i=0;i<size;i++){
		total_proTH1+=(TH0array[i].proTH1/100);
		total_proTH2+=(TH0array[i].proTH2/100);
	}
	oxyDeficit=(xDim*yDim)-system_oxy;
//	if(step==1745){cout<<"ProcID="<<procID<<" oxy="<<oxyDeficit<<" Infect="<<totalInfection<<"\n";cout.flush();}
	
}

void readParams(int procID, int nProc,int p1[], int* np1, float p2[], int p3[], int p4[], 
	int p5[], int p6[]){
	
	ifstream paramFile1,paramFile2;
	int intParams[parameterInput][3]; //NIR,IS,NRI
	float floatParams[parameterInput]; //healParam
	int i,j,ic,index;
	int inj_number,seed,newID,counter;
	
	paramFile1.open("NIP.dat");
	paramFile2.open("NFP.dat");
	i=0;
	j=0;
	index=0;
	newID=procID+1;
	counter=1;
	while(!paramFile1.eof()){
		paramFile1>>intParams[i][j];
		j++;
		if(j==3){
			j=0;
			i++;}
		if(i==parameterInput){break;}
	}
	i=0;
	while(!paramFile2.eof()){
		paramFile2>>floatParams[i];
		i++;
		if(i==parameterInput){break;}
	}
	paramFile1.close();
	paramFile2.close();
					
	for(inj_number=inj_numberStart;inj_number<=inj_numberMax;inj_number+=inj_numberIncrement){
	for(seed=seedStart;seed<=seedMax;seed+=seedIncrement){
	for(ic=0;ic<parameterInput;ic++){
		if(counter==newID){
			p1[index]=inj_number;
			(*np1)++;
			p2[index]=floatParams[ic];
			p3[index]=intParams[ic][1];
			p4[index]=intParams[ic][2];
			p5[index]=intParams[ic][0];
			p6[index]=seed;
			index++;
		}
		counter=counter+1;
		if(counter>nProc){counter=1;}
		
	}}}
	cout<<"ProcID="<<procID<<" Parameter Distribution Complete\n";
}

void pickParams(int p1[], int* np1, float p2[], int p3[], int p4[], int p5[], int p6[], 
	int procID, int maxProcs){

	int counter,newID,index;
	
	newID=procID+1;
	counter=1;
	index=0;
	
	float oxyHeal;
      int inj_number,numRecurInj,infectSpread,numInfectRepeat,seed;

	for(inj_number=inj_numberStart;inj_number<=inj_numberMax;inj_number+=inj_numberIncrement){
	for(oxyHeal=oxyHealStart;oxyHeal<=oxyHealMax;oxyHeal+=oxyHealIncrement){
	for(infectSpread=infectSpreadStart;infectSpread<=infectSpreadMax;infectSpread+=infectSpreadIncrement){
	for(numRecurInj=numRecurInjStart;numRecurInj<=numRecurInjMax;numRecurInj+=numRecurInjIncrement){
	for(numInfectRepeat=numInfectRepeatStart;numInfectRepeat<=numInfectRepeatMax;numInfectRepeat+=numInfectRepeatIncrement){
	for(seed=seedStart;seed<=seedMax;seed+=seedIncrement){
		if(counter==newID){
			p1[index]=inj_number;
			(*np1)++;
			p2[index]=oxyHeal;
			p3[index]=infectSpread;
			p4[index]=numRecurInj;
			p5[index]=numInfectRepeat;
			p6[index]=seed;
			index++;
		}
		counter=counter+1;
		if(counter>maxProcs){counter=1;}
	}}}}}}
}

void adjustOrientation(int* orientation, int leftOrRight){
//if leftOrRight=-1, adjust orientation left, if =1, adjust orientation right
	int tempOrient;
	tempOrient=*orientation;
	tempOrient+=leftOrRight;
	if(tempOrient>7){tempOrient=0;}
	if(tempOrient<0){tempOrient=7;}
	*orientation=tempOrient;
}

void move(int orient, int* x, int* y){
	int oldx,newx,oldy,newy;
	
	oldx=*x;
	oldy=*y;
	newx=oldx;
	newy=oldy;

	if(orient==0){ //Move North
		newy=oldy+1;
		if(newy>=yDim){newy=0;}
	}
	
	if(orient==1){ //Move Northeast
		newy=oldy+1;
		newx=oldx+1;
		if(newy>=yDim){newy=0;}
		if(newx>=yDim){newx=0;}
	}
	
	if(orient==2){ //Move East
		newx=oldx+1;
		if(newx>=xDim){newx=0;}
	}
	
	if(orient==3){ //Move Southeast
		newy=oldy-1;
		newx=oldx+1;
		if(newx>=xDim){newx=0;}
		if(newy<0){newy=yDim-1;}
	}
	
	if(orient==4){ //Move South
		newy=oldy-1;
		if(newy<0){newy=yDim-1;}
	}
	
	if(orient==5){ //Move Southwest
		newy=oldy-1;
		newx=oldx-1;
		if(newy<0){newy=yDim-1;}
		if(newx<0){newx=xDim-1;}
	}
	
	if(orient==6){ //Move West
		newx=oldx-1;
		if(newx<0){newx=xDim-1;}
	}
	
	if(orient==7){ //Move Northwest
		newx=oldx-1;
		newy=oldy+1;
		if(newx<0){newx=xDim-1;}
		if(newy>=yDim){newy=0;}
	}

	if(cellGrid[newx][newy]<cellCapacity){
		*x=newx;
		*y=newy;
		cellGrid[oldx][oldy]--;
		cellGrid[newx][newy]++;
	}	
}

void getAhead(int orient, int x, int y, int *xl, int *xm, int *xr, int *yl, int *ym, int *yr){
	int txl,txm,txr,tyl,tym,tyr;

	if(orient==0){
		txl=x-1;
		txm=x;
		txr=x+1;
		tyl=y+1;
		tym=y+1;
		tyr=y+1;
	}
	
	if(orient==1){
		txl=x;
		txm=x+1;
		txr=x+1;
		tyl=y+1;
		tym=y+1;
		tyr=y;
	}
	
	if(orient==2){
		txl=x+1;
		txm=x+1;
		txr=x+1;
		tyl=y+1;
		tym=y;
		tyr=y-1;
	}
	
	if(orient==3){
		txl=x+1;
		txm=x+1;
		txr=x;
		tyl=y;
		tym=y-1;
		tyr=y-1;
	}
	
	if(orient==4){
		txl=x+1;
		txm=x;
		txr=x-1;	
		tyl=y-1;
		tym=y-1;
		tyr=y-1;
	}
	
	if(orient==5){
		txl=x;
		txm=x-1;
		txr=x-1;
		tyl=y-1;
		tym=y-1;
		tyr=y;
	}
	
	if(orient==6){
		txl=x-1;
		txm=x-1;
		txr=x-1;
		tyl=y-1;
		tym=y;
		tyr=y+1;
	}
	
	if(orient==7){
		txl=x-1;
		txm=x-1;
		txr=x;
		tyl=y;
		tym=y+1;
		tyr=y+1;
	}
	
	if(txl<0){txl=xDim-1;}
	if(txl>=xDim){txl=0;}
	if(txm<0){txm=xDim-1;}
	if(txm>=xDim){txm=0;}
	if(txr<0){txr=xDim-1;}
	if(txr>=xDim){txr=0;}
	if(tyl<0){tyl=yDim-1;}
	if(tyl>=yDim){tyl=0;}
	if(tym<0){tym=yDim-1;}
	if(tym>=yDim){tym=0;}
	if(tyr<0){tyr=yDim-1;}
	if(tyr>=yDim){tyr=0;}

	*xl=txl;
	*xm=txm;
	*xr=txr;
	*yl=tyl;
	*ym=tym;
	*yr=tyr;
}

void wiggle(int* orientation){ //Should always be followed by move() to match NetLogo
	int dir,tempOrient;
	tempOrient=*orientation;
	dir=distribution3(generator);
	if(dir==0){tempOrient--;}
	if(dir==2){tempOrient++;}
	if(tempOrient>7){tempOrient=0;}
	if(tempOrient<0){tempOrient=7;}
	*orientation=tempOrient;
}

void injure_infectionFRD(int inj_number){ //Fixed Radius Disk
	int i,x,y,rad,size,id;
	size=ecArray.size();
	rad=inj_number;
	for(i=0;i<size;i++){
		x=ecArray[i].xLoc-xDim/2;
		y=ecArray[i].yLoc-yDim/2;
		if((pow(x,2)+pow(y,2))<=(pow(rad,2))){
			id=ecArray[i].yLoc*xDim+ecArray[i].xLoc;
			ecArray[id].infection=100;
		}
	}
}

void recur_injury(){
	int x,y,id;
	x=distribution100(generator);
	y=distribution100(generator);
	id=y*xDim+x;
	ecArray[id].infection=100;
}

void EC::getNeighbors(){
	int nw,n,ne,e,se,s,sw,w;

	nw=id+xDim-1;
	n=id+xDim;
	ne=id+xDim+1;
	e=id+1;
	se=id-xDim+1;
	s=id-xDim;
	sw=id-xDim-1;
	w=id-1;
	
	if(id%xDim==0){ //object lies on west grid border
		nw=id+2*xDim-1;
		w=id+xDim-1;
		sw=id-1;
	}
	
	if(id%xDim==xDim-1){ //object lies on east grid border
		ne=id+1;
		e=id-(xDim-1);
		se=id-2*xDim+1;
	}
	
	if(id+xDim>=xDim*yDim){ //object lies on northern border
		n=(id+xDim)%(xDim*yDim);
		ne=(id+xDim)%(xDim*yDim)+1;
		nw=(id+xDim)%(xDim*yDim)-1;
		if(id%xDim==0){
			nw=(id+xDim)%(xDim*yDim)+xDim-1;
		}
		if(id%xDim==xDim-1){
			ne=(id+xDim)%(xDim*yDim)-(xDim-1);
		}	
	}
	
	if(id<xDim){ //object lies on southern border
		s=(xDim*yDim)-xDim+id;
		se=(xDim*yDim)-xDim+id+1;
		sw=(xDim*yDim)-xDim+id-1;
		if(id%xDim==0){
			sw=(xDim*yDim)-1;
		}
		if(id%xDim==xDim-1){
			se=(xDim*yDim)-xDim;
		}	
	}
	
	neighbors[0]=nw;
	neighbors[1]=n;
	neighbors[2]=ne;
	neighbors[3]=e;
	neighbors[4]=se;
	neighbors[5]=s;
	neighbors[6]=sw;
	neighbors[7]=w;
}

void initialize(){
	int i,j,k,xTemp,yTemp;
	
	ecArray.clear();      
	pmnArray.clear();
	monoArray.clear();
	TH0array.clear();
	TH1array.clear();
	TH2array.clear();
	pmn_marrowArray.clear();
	mono_marrowArray.clear();
	TH0_germArray.clear();
	TH1_germArray.clear();
	TH2_germArray.clear();

	system_oxy=10201;
	oxyDeficit=0;
	totalInfection=0;
	total_TNF=0;
	total_sTNFr=0;
	total_IL10=0;
	total_GCSF=0;
	total_proTH1=0;
	total_proTH2=0;
	
	PAFmult=1;
	TNFmult=1;
	sTNFrmult=1;
	IL1mult=1;
	sIL1rmult=1;
	IL1ramult=1;
	IFNgmult=1;
	IL4mult=1;
	IL8mult=1;
	IL10mult=1;
	IL12mult=1;
	GCSFmult=1;

	mGCSF=0;
	mPAF=0;
	mTNF=0;
	mSTNFR=0;
	mIL1=0;
	mSIL1R=0;
	mIL1RA=0;
	mIFNg=0;
	mIL4=0;
	mIL8=0;
	mIL10=0;
	mIL12=0;
	
	for(i=0;i<xDim;i++){
		for(j=0;j<yDim;j++){
			cellGrid[i][j]=0;}}
	
	k=0; //initialization
	for(j=0;j<yDim;j++){       //Initialize EC grid
		for(i=0;i<xDim;i++){
			ecArray.push_back(EC(i,j,k));
			ecArray[k].getNeighbors();
			k++;
		}
	}

	for(i=0;i<500;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
		pmnArray.push_back(pmn(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
	
	for(i=0;i<50;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		monoArray.push_back(mono(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
	
	for(i=0;i<50;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		TH1array.push_back(TH1(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
	
	for(i=0;i<50;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		TH2array.push_back(TH2(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		pmn_marrowArray.push_back(pmn_marrow(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		mono_marrowArray.push_back(mono_marrow(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
		
	for(i=0;i<100;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		TH0_germArray.push_back(TH0_germ(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		TH1_germArray.push_back(TH1_germ(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
	
	for(i=0;i<100;i++){
		xTemp=distribution100(generator);
		yTemp=distribution100(generator);
//		cout<<xTemp<<" "<<yTemp<<"\n";
		TH2_germArray.push_back(TH2_germ(xTemp,yTemp));
		cellGrid[xTemp][yTemp]++;
	}
}

void spreadParams(int p1[], float p2[], int p3[], int p4[], int p5[], int p6[], int* np1){
			p1[0]=inj_numberStart;
			p2[0]=oxyHealStart;
			p3[0]=infectSpreadStart;
			p4[0]=numRecurInjStart;
			p5[0]=numInfectRepeatStart;
			p6[0]=seedStart;
			(*np1)++;
}

void applyIntervention(int flag){
	PAFmult=myIntervention[flag*12+0];
	TNFmult=myIntervention[flag*12+1];
	sTNFrmult=myIntervention[flag*12+2];
	IL1mult=myIntervention[flag*12+3];
	sIL1rmult=myIntervention[flag*12+4];
	IL1ramult=myIntervention[flag*12+5];
	IFNgmult=myIntervention[flag*12+6];
	IL4mult=myIntervention[flag*12+7];
	IL8mult=myIntervention[flag*12+8];
	IL10mult=myIntervention[flag*12+9];
	IL12mult=myIntervention[flag*12+10];
	GCSFmult=myIntervention[flag*12+11];		
}

void applyTestIntervention(int flag){
//	cout<<"Intervention Applied "<<flag<<"\n";
	PAFmult=testInt[flag*12+0];
	TNFmult=testInt[flag*12+1];
	sTNFrmult=testInt[flag*12+2];
	IL1mult=testInt[flag*12+3];
	sIL1rmult=testInt[flag*12+4];
	IL1ramult=testInt[flag*12+5];
	IFNgmult=testInt[flag*12+6];
	IL4mult=testInt[flag*12+7];
	IL8mult=testInt[flag*12+8];
	IL10mult=testInt[flag*12+9];
	IL12mult=testInt[flag*12+10];
	GCSFmult=testInt[flag*12+11];
}

void clearIntervention(){
	PAFmult=1;
	TNFmult=1;
	sTNFrmult=1;
	IL1mult=1;
	sIL1rmult=1;
	IL1ramult=1;
	IFNgmult=1;
	IL4mult=1;
	IL8mult=1;
	IL10mult=1;
	IL12mult=1;
	GCSFmult=1;

}

void readIntervention(){
	ifstream interventionFile;
	int i=0;
	interventionFile.open(intFile);
	if(testIntervention==1&&GA==0){
		while(!interventionFile.eof()){
			interventionFile>>myIntervention[i];
			i++;
		}
	}
	if(testIntervention==1&&GA==1){
		while(!interventionFile.eof()){
			interventionFile>>testInt[i];
			i++;
		}
	}

}
