#include <vector>
#include "agents.h"
#include <stdlib.h>
#include <random>
using namespace std;

extern mt19937 generator;
extern uniform_int_distribution<int> distribution8,distribution50;

EC::EC(){}

EC::EC(int x, int y, int iid){
	xLoc=x;
	yLoc=y;
	id=iid;
	oxy=100;
	ec_activation=0;
	ec_roll=0;
	ec_stick=0;
	ec_migrate=0;
	endotoxin=0;
	infection=0;
	cytotox=0;
	PAF=0;
	TNF=0;
	sTNFr=0;
	IL1=0;
	sIL1r=0;
	IL1ra=0;
	IFNg=0;
	IL4=0;
	IL8=0;
	IL10=0;
	IL12=0;
	GCSF=0;

}

pmn::pmn(int x, int y){
	xLoc=x;
	yLoc=y;
	pmn_age=distribution50(generator);
	wbc_roll=1;
	wbc_stick=0;
	wbc_migrate=0;
	pmn_pcd=10;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<pmn_age<<" "<<orientation<<"\n";
}

pmn::pmn(int x, int y, int age){
	xLoc=x;
	yLoc=y;
	pmn_age=age;
	wbc_roll=1;
	wbc_stick=0;
	wbc_migrate=0;
	pmn_pcd=10;
	orientation=distribution8(generator); //n,e,s,w
}

mono::mono(int x, int y){
	xLoc=x;
	yLoc=y;
	mono_age=distribution50(generator);
	TNFr=0;
	IL_1r=0;
	activation=0;
	wbc_roll=1;
	wbc_stick=0;
	wbc_migrate=0;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<mono_age<<" "<<orientation<<"\n";
}

mono::mono(int x, int y, int age, int iTNFr, int iIL1r){
	xLoc=x;
	yLoc=y;
	mono_age=age;
	TNFr=iTNFr;
	IL_1r=iIL1r;
	activation=0;
	orientation=distribution8(generator); //n,e,s,w
	wbc_roll=1;
	wbc_stick=0;
	wbc_migrate=0;
}

TH0::TH0(int x, int y){
	xLoc=x;
	yLoc=y;
	TH0_age=distribution50(generator);
	proTH1=0;
	proTH2=0;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<TH0_age<<" "<<orientation<<"\n";
}

TH0::TH0(int x, int y, int age){
	xLoc=x;
	yLoc=y;
	TH0_age=age;
	proTH1=0;
	proTH2=0;
	orientation=distribution8(generator); //n,e,s,w
}

TH1::TH1(int x, int y){
	xLoc=x;
	yLoc=y;
	TH1_age=distribution50(generator);
	orientation=distribution8(generator); //n,e,s,w
//	cout<<TH1_age<<" "<<orientation<<"\n";
}

TH1::TH1(int x, int y, int age){
	xLoc=x;
	yLoc=y;
	TH1_age=age;
	orientation=distribution8(generator); //n,e,s,w
}

TH2::TH2(int x, int y){
	xLoc=x;
	yLoc=y;
	TH2_age=distribution50(generator);
	orientation=distribution8(generator); //n,e,s,w
//	cout<<TH2_age<<" "<<orientation<<"\n";
}

TH2::TH2(int x, int y, int age){
	xLoc=x;
	yLoc=y;
	TH2_age=age;
	orientation=distribution8(generator); //n,e,s,w
}

pmn_marrow::pmn_marrow(int x, int y){
	xLoc=x;
	yLoc=y;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<orientation<<"\n";
}

mono_marrow::mono_marrow(int x, int y){
	xLoc=x;
	yLoc=y;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<orientation<<"\n";
}

TH0_germ::TH0_germ(int x, int y){
	xLoc=x;
	yLoc=y;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<orientation<<"\n";
}

TH1_germ::TH1_germ(int x, int y){
	xLoc=x;
	yLoc=y;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<orientation<<"\n";
}

TH2_germ::TH2_germ(int x, int y){
	xLoc=x;
	yLoc=y;
	orientation=distribution8(generator); //n,e,s,w
//	cout<<orientation<<"\n";
}


