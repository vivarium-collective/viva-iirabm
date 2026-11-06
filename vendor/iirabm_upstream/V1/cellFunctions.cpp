#include <stdlib.h>
#include <vector>
#include <iostream>
#include "agents.h"
#include <algorithm>
#include <cmath>
#include "Parameters.h"
using namespace std;

extern int cellGrid[101][101],procID;
extern float total_GCSF;
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
	
extern mt19937 generator;
extern uniform_int_distribution<int> distribution10k,distribution100,distribution10,distribution8;

void getAhead(int orient, int x, int y, int *xl, int *xm, int *xr, int *yl, int *ym, int *yr);
extern void wiggle(int* orientation);
extern void adjustOrientation(int* orientation, int leftOrRight);
extern void move(int orient, int* x, int* y);
extern void injure_infectionFRD(int inj_number);	

void heal(int index);
void applyAntibiotics();

void EC::inj_function(int infectSpread, int numInfectRepeat){
	int temp,i;
	
	oxy=max(float(0),oxy-infection);
	endotoxin+=(infection/10);
	
	for(i=1;i<=numInfectRepeat;i++){
		if(infection>=100){
			temp=distribution8(generator);
			ecArray[neighbors[temp]].infection+=infectSpread;
			infection=100;
	}}
	if(infection>0){
		infection=max(float(0),infection-cytotox+float(0.1));}
}

void EC::activate(){
	ec_roll++;
	ec_stick++;
	if(PAFmult<=1){
		PAF=(PAF+1)*PAFmult;
	} else{
		PAF=PAF+1+PAFmult-1;
	}
	if(IL8mult<=1){	
		IL8=(IL8+1)*IL8mult;
	} else{
		IL8=IL8+1+IL8mult-1;
	}
}

void EC::ECfunction(float oxyHeal){
	if((endotoxin>=1)||(oxy<60)){
		ec_activation=1;}
	if(ec_activation==1){
		activate();}
	patch_inj_spread(oxyHeal);
}

void EC::patch_inj_spread(float oxyHeal){
	int i;
	oxy=oxy-cytotox;
	
	if(oxy>=60){
		oxy=min(float(100),oxy+oxyHeal);}
	
	if((oxy<60)&&(oxy>30)){ //ischemia
		ec_roll++;
		oxy-=0.5;
		if(PAFmult<=1){
			PAF=(PAF+1)*PAFmult;
		} else{
			PAF=PAF+1+PAFmult-1;
		}
		for(i=0;i<8;i++){
			ecArray[neighbors[i]].oxy-=0.05;}
	}
		if(oxy<=30){ //infarction
			ec_stick++;
			oxy-=2;
			if(PAFmult<=1){
				PAF=(PAF+1)*PAFmult;
			} else{
				PAF=PAF+1+PAFmult-1;
			}
		for(i=0;i<8;i++){
			ecArray[neighbors[i]].oxy-=0.25;
			if(ecArray[neighbors[i]].oxy<0){ecArray[neighbors[i]].oxy=0;}
		}
	}
	if(oxy<0){oxy=0;}
}

void pmn::pmn_function(int pmnID){
	int roll,x,y,id;
	float TNF,PAF,IL1,IL10,GCSF;
	x=xLoc;
	y=yLoc;
	id=y*xDim+x;
	roll=ecArray[id].ec_roll;
	TNF=ecArray[id].TNF;
	PAF=ecArray[id].PAF;
	IL1=ecArray[id].IL1;
	IL10=ecArray[id].IL10;
	GCSF=ecArray[id].GCSF;
	if(wbc_migrate>0){
		pmn_burst(pmnID);}
	else{
		if((roll>3)&&(wbc_roll==1)){
			pmn_sniff();}
		else{
			pmn_sniff();
			pmn_sniff();}
		x=xLoc;
		y=yLoc;
		id=y*xDim+x;
		TNF=ecArray[id].TNF;
		PAF=ecArray[id].PAF;
		IL1=ecArray[id].IL1;
		IL10=ecArray[id].IL10;
		GCSF=ecArray[id].GCSF;
		if(TNF+PAF>1){
			wbc_stick=IL1;
			if(IL1ramult<=1){
				ecArray[id].IL1ra=(ecArray[id].IL1ra+1)*IL1ramult;
			} else{
				ecArray[id].IL1ra=ecArray[id].IL1ra+1+IL1ramult-1;
			}
		}
		if((wbc_stick>=1)&&(ecArray[id].ec_stick>=100)){
			wbc_migrate=max(float(0),(TNF+IL1+GCSF-IL10));}
		pmn_age--;
		if(pmn_age<0){
			x=xLoc;
			y=yLoc;
			pmnArray.erase(pmnArray.begin()+pmnID);
			cellGrid[x][y]--;
			if(cellGrid[x][y]<0){cout<<"Cell Error PMNFUNC\n";
			return;
			}}
	}
}

void pmn::pmn_burst(int pmnID){
	int x,y,id;
	
	x=xLoc;
	y=yLoc;
	id=y*xDim+x;
	ecArray[id].cytotox=max(float(10),ecArray[id].TNF);
	ecArray[id].oxy=100;
	ecArray[id].ec_roll=0;
	ecArray[id].ec_stick=0;
	ecArray[id].ec_migrate=0;
	if(TNFmult<=1){
		ecArray[id].TNF=(ecArray[id].TNF+1)*TNFmult;
	} else{
		ecArray[id].TNF=ecArray[id].TNF+1+TNFmult-1;
	}
	if(IL1mult<=1){
		ecArray[id].IL1=(ecArray[id].IL1+1)*IL1mult;
	} else{
		ecArray[id].IL1=ecArray[id].IL1+1+IL1mult-1;
	}
	pmn_age=pmn_pcd;
	pmn_pcd=pmn_pcd-1+max(float(0),(ecArray[id].TNF+ecArray[id].IFNg+ecArray[id].GCSF-
		ecArray[id].IL10)/100);
	if(pmn_age<0){
		x=xLoc;
		y=yLoc;
		pmnArray.erase(pmnArray.begin()+pmnID);
		cellGrid[x][y]--;
	}
}

void pmn::pmn_sniff(){
	int x,y,idm,idr,idl,xl,xm,xr,yl,ym,yr,flag;
	float pmnahead,pmnright,pmnleft;
	x=xLoc;
	y=yLoc;

	flag=0;  //Flag=-1 for left, 0 for middle, and 1 for right

	getAhead(orientation, x, y, &xl, &xm, &xr, &yl, &ym, &yr);
	
	idr=yr*xDim+xr;
	idm=ym*xDim+xm;
	idl=yl*xDim+xl;
	pmnright=ecArray[idr].IL8;
	pmnahead=ecArray[idm].IL8;
	pmnleft=ecArray[idl].IL8;

	if((pmnright>=pmnahead)&&(pmnright>=pmnleft)){flag=1;}
	else if(pmnleft>=pmnahead){flag-=1;}
	adjustOrientation(&orientation,flag);
	
	if(flag==-1){
		if(cellGrid[xl][yl]<cellCapacity){
 			xLoc=xl;
 			yLoc=yl;
 			cellGrid[xl][yl]++;
 			cellGrid[x][y]--;}
	}
	if(flag==0){
		if(cellGrid[xm][ym]<cellCapacity){
 			xLoc=xm;
 			yLoc=ym;
 			cellGrid[xm][ym]++;
 			cellGrid[x][y]--;}
	}
	if(flag==1){
		if(cellGrid[xr][yr]<cellCapacity){
 			xLoc=xr;
 			yLoc=yr;
 			cellGrid[xr][yr]++;
 			cellGrid[x][y]--;}
	}
}

void mono::mono_function(int index){
	int x,y,id;
	
	x=xLoc;
	y=yLoc;
	id=y*xDim+x;
	
	if(ecArray[id].sTNFr<=100){
		TNFr=min(float(100),(ecArray[id].TNF+ecArray[id].sTNFr));}
	else{
		TNFr=min(float(100),max(float(0),ecArray[id].TNF-ecArray[id].sTNFr));}
	IL_1r=min(float(100),max(float(0),ecArray[id].IL1-ecArray[id].IL1ra-ecArray[id].sIL1r));
	
	if(IL1ramult<=1){
		ecArray[id].IL1ra=(ecArray[id].IL1ra+ecArray[id].IL1/2)*IL1ramult;
	} else {
		ecArray[id].IL1ra=(ecArray[id].IL1ra+ecArray[id].IL1/2)+IL1ramult-1;
	}
	if(sTNFrmult<=1){
		ecArray[id].sTNFr=(ecArray[id].sTNFr+TNFr/2)*sTNFrmult;
	} else {
		ecArray[id].sTNFr=(ecArray[id].sTNFr+TNFr/2)+sTNFrmult-1;
	}
	if(sIL1rmult<=1){
		ecArray[id].sIL1r=(ecArray[id].sIL1r+IL_1r/2)*sIL1rmult;
	} else {
		ecArray[id].sIL1r=(ecArray[id].sIL1r+IL_1r/2)+sIL1rmult-1;
	}
	activation=ecArray[id].endotoxin+ecArray[id].PAF+ecArray[id].IFNg-ecArray[id].IL10;
	if(activation>0){
		if(GCSFmult<=1){
			ecArray[id].GCSF=(ecArray[id].GCSF+ecArray[id].endotoxin+ecArray[id].PAF+ecArray[id].TNF+
			ecArray[id].IFNg)*GCSFmult;
		} else {
			ecArray[id].GCSF=(ecArray[id].GCSF+ecArray[id].endotoxin+ecArray[id].PAF+ecArray[id].TNF+
			ecArray[id].IFNg)+GCSFmult-1;
		}
		if(IL8mult<=1){
			ecArray[id].IL8=(ecArray[id].IL8+ecArray[id].TNF+ecArray[id].IL1)*IL8mult;
		} else {
			ecArray[id].IL8=(ecArray[id].IL8+ecArray[id].TNF+ecArray[id].IL1)+IL8mult-1;
		}
		if(IL12mult<=1){
			ecArray[id].IL12=(ecArray[id].IL12+ecArray[id].TNF+ecArray[id].IL1)*IL12mult;
		} else {
			ecArray[id].IL12=(ecArray[id].IL12+ecArray[id].TNF+ecArray[id].IL1)+IL12mult;
		}
		if(IL10mult<=1){
			ecArray[id].IL10=(ecArray[id].IL10+ecArray[id].TNF+ecArray[id].IL1)*IL10mult;
		} else {
			ecArray[id].IL10=(ecArray[id].IL10+ecArray[id].TNF+ecArray[id].IL1)+IL10mult-1;
		}
		if(IL1mult<=1){
			ecArray[id].IL1=(ecArray[id].IL1+ecArray[id].endotoxin+ecArray[id].PAF+IL_1r+
			ecArray[id].TNF)*IL1mult;
		} else { 
			ecArray[id].IL1=(ecArray[id].IL1+ecArray[id].endotoxin+ecArray[id].PAF+IL_1r+
			ecArray[id].TNF)+IL1mult-1;
		}
		if(TNFmult<=1){
			ecArray[id].TNF=(ecArray[id].TNF+ecArray[id].endotoxin+ecArray[id].PAF+TNFr+
			ecArray[id].IFNg)*TNFmult;	
		} else {
			ecArray[id].TNF=(ecArray[id].TNF+ecArray[id].endotoxin+ecArray[id].PAF+TNFr+
			ecArray[id].IFNg)+TNFmult-1;
		}
		if((wbc_stick==1)&&(ecArray[id].ec_stick>=100)){
			wbc_migrate=1;}
		if(wbc_roll==1){
			wbc_stick=1;}
		wbc_roll=1;
	}
	if(activation<0){
		if(IL10mult<=1){
			ecArray[id].IL10=(ecArray[id].IL10+ecArray[id].TNF+ecArray[id].IL1)*IL10mult;
		} else {
			ecArray[id].IL10=(ecArray[id].IL10+ecArray[id].TNF+ecArray[id].IL1)+IL10mult-1;
		}
	}
	if(wbc_migrate==1){
		heal(id);}
	if(wbc_roll==1){
		mono_sniff();}
	else{
		mono_sniff();
		mono_sniff();}
	mono_age--;
	if(mono_age<0){
		x=xLoc;
		y=yLoc;
		monoArray.erase(monoArray.begin()+index);
		cellGrid[x][y]--;
	}
	if(activation>20){
		activation=20;}
}

void mono::mono_sniff(){
	int x,y,idm,idr,idl,xl,xm,xr,yl,ym,yr,flag;
	float pafahead,pafright,pafleft;
	x=xLoc;
	y=yLoc;
	
	flag=0;
	getAhead(orientation, x, y, &xl, &xm, &xr, &yl, &ym, &yr);
	idr=yr*xDim+xr;
	idm=ym*xDim+xm;
	idl=yl*xDim+xl;
	pafright=ecArray[idr].PAF;
	pafahead=ecArray[idm].PAF;
	pafleft=ecArray[idl].PAF;

	if((pafright>=pafahead)&&(pafright>=pafleft)){flag=1;}
	else if(pafleft>=pafahead){flag-=1;}
	adjustOrientation(&orientation,flag);
	
	if(flag==-1){
		if(cellGrid[xl][yl]<cellCapacity){
 			xLoc=xl;
 			yLoc=yl;
 			cellGrid[xl][yl]++;
 			cellGrid[x][y]--;}
	}
	
	if(flag==0){
		if(cellGrid[xm][ym]<cellCapacity){
 			xLoc=xm;
 			yLoc=ym;
 			cellGrid[xm][ym]++;
 			cellGrid[x][y]--;}
	}
	
	if(flag==1){
		if(cellGrid[xr][yr]<cellCapacity){
 			xLoc=xr;
 			yLoc=yr;
 			cellGrid[xr][yr]++;
 			cellGrid[x][y]--;}
	}
}

void TH0::TH0function(int index){
	int id;
	id=yLoc*xDim+xLoc;
	if(ecArray[id].IL12+ecArray[id].IL4>0){
		proTH1=(ecArray[id].IL12+ecArray[id].IFNg)*100;
		proTH2=(ecArray[id].IL10+ecArray[id].IL4)*100;
		if((proTH1>0)&&(proTH2>0)){
			rTH1=distribution10k(generator)%int(ceil(proTH1));
			rTH2=distribution10k(generator)%int(ceil(proTH2));
			if(rTH1>rTH2){activation++;}
			if(rTH1<rTH2){activation--;}
		}
		if(proTH1==0){activation--;}
		if(proTH2==0){activation++;}
	}
	wiggle(&orientation);	
	move(orientation, &xLoc, &yLoc);
	TH0_age--;
	if(TH0_age<0){
		cellGrid[xLoc][yLoc]--;
		TH0array.erase(TH0array.begin()+index);
		return;}
	if(activation>=10){
		TH1array.push_back(TH1(xLoc,yLoc,TH0_age));
		TH0array.erase(TH0array.begin()+index);
	}
	if(activation<=-10){
		TH2array.push_back(TH2(xLoc,yLoc,TH0_age));
		TH0array.erase(TH0array.begin()+index);
	}
}

void TH1::TH1function(int index){
	int id;

	id=yLoc*xDim+xLoc;
	if(ecArray[id].IL12>0){
		if(IFNgmult<=1){
			ecArray[id].IFNg=(2*ecArray[id].IFNg+ecArray[id].IL12+ecArray[id].TNF+
			ecArray[id].IL1)*IFNgmult;
		} else {
			ecArray[id].IFNg=(2*ecArray[id].IFNg+ecArray[id].IL12+ecArray[id].TNF+
			ecArray[id].IL1)+IFNgmult-1;
		}
	}
	wiggle(&orientation);
	move(orientation, &xLoc, &yLoc);
	TH1_age--;
	if(TH1_age<0){
		cellGrid[xLoc][yLoc]--;
		TH1array.erase(TH1array.begin()+index);
	}
}

void TH2::TH2function(int index){
	int id;
	
	id=yLoc*xDim+xLoc;
	
	if(ecArray[id].IL10>0){
		if(IL4mult<=1){
			ecArray[id].IL4=(ecArray[id].IL4+ecArray[id].IL10)*IL4mult;
		} else {
			ecArray[id].IL4=(ecArray[id].IL4+ecArray[id].IL10)+IL4mult-1;
		}
		if(IL10mult<=1){
			ecArray[id].IL10=ecArray[id].IL10*2*IL10mult;
		} else {
			ecArray[id].IL10=ecArray[id].IL10*2+IL10mult-1;
		}
	}
	wiggle(&orientation);
	move(orientation, &xLoc, &yLoc);
	TH2_age--;
	if(TH2_age<0){
		cellGrid[xLoc][yLoc]--;
		TH2array.erase(TH2array.begin()+index);
		
	}
}

void pmn_marrow::pmn_marrow_function(){
	int x,y,temp,n,i;
	if(total_GCSF>3000){total_GCSF=3000;}
	n=int(1+total_GCSF/100);
	
	for(i=0;i<n;i++){
		x=distribution100(generator);
		y=distribution100(generator);
		temp=distribution10(generator);
		if(temp<1){
			pmnArray.push_back(pmn(x,y,50));
			cellGrid[x][y]++;
		}
	}
}

void mono_marrow::mono_marrow_function(){
	int x,y,temp;
	x=distribution100(generator);
	y=distribution100(generator);
	temp=distribution100(generator);
	if(temp<1){
		monoArray.push_back(mono(x,y,50,1000,1000));
		cellGrid[x][y]++;
	}
}

void TH0_germ::TH0_germ_function(){
	int x,y,temp;
	x=xLoc;
	y=yLoc;
	temp=distribution100(generator);
	if(temp<1){
		TH0array.push_back(TH0(x,y,100));
		cellGrid[x][y]++;
	}
}

void TH1_germ::TH1_germ_function(){
	int x,y,temp;
	x=xLoc;
	y=yLoc;
	temp=distribution100(generator);
	if(temp<1){
		TH1array.push_back(TH1(x,y,100));
		cellGrid[x][y]++;
	}
}

void TH2_germ::TH2_germ_function(){
	int x,y,temp;
	x=xLoc;
	y=yLoc;
	temp=distribution100(generator);
	if(temp<1){
		TH2array.push_back(TH2(x,y,100));
		cellGrid[x][y]++;
	}
}

void heal(int index){
	ecArray[index].oxy=100;
	ecArray[index].ec_roll=0;
	ecArray[index].ec_stick=0;
	ecArray[index].ec_migrate=0;
	ecArray[index].infection=0;
	ecArray[index].ec_activation=0;
}

void applyAntibiotics(){
	int i,length;
	length=ecArray.size();
	
	for(i=0;i<length;i++){
		ecArray[i].infection=ecArray[i].infection*antibioticMultiplier;	
	}
}

