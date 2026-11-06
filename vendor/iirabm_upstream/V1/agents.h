#ifndef AGENTS_H
#define AGENTS_H
#include<vector>
#include <stdlib.h>
#include <iostream>
using namespace std;

struct EC{
   public:
	int xLoc;
	int yLoc;
	int id;
	bool alive;
	float oxy;
	int ec_activation;
	int ec_roll;         //rolling
	int ec_stick;        //sticking
	int ec_migrate;      //migration
	float cytotox;         //o2rads and enzymes
	float infection;       //infectious vector
	float endotoxin;       
	float PAF;             
	float TNF;
	float sTNFr;
	float IL1;
	float sIL1r;
	float IL1ra;
	float IFNg;
	float IL4;
	float IL6;
	float IL8;
	float IL10;
	float IL12;
	float GCSF; 

	float tempSignal;
	int neighbors[8];
	EC(int x, int y, int id);
	EC();
	void ECfunction(float oxyHeal);
	void inj_function(int infectSpread, int numInfectRepeat);
	void getNeighbors();
	void diffuse();
	void activate();
	void patch_inj_spread(float oxyHeal);
};

struct pmn{
   public:
	int xLoc;
	int yLoc;
	int wbc_roll;        //selectins
	int wbc_stick;       //integrens
	int wbc_migrate;     //diapedesis
	int pmn_age;         //life span
	float pmn_pcd;         //programmed cell death
	float TNFr;
	float IL_1r;
	float activation;
	int orientation;
	pmn(int x, int y);
	pmn(int x, int y, int age);
	void pmn_function(int pmnID);
	void pmn_burst(int pmnID);
	void pmn_sniff();
};

struct pmn_marrow{
   public:
	int xLoc;
	int yLoc;
	float wbc_roll;        //selectins
	float wbc_stick;       //integrens
	float wbc_migrate;     //diapedesis
	int pmn_age;         //life span
	float pmn_pcd;         //programmed cell death
	int mono_age;        //life span
	float TNFr;
	float IL_1r;
	float activation;
	pmn_marrow(int x, int y);
	void pmn_marrow_function();
	int orientation;
};

struct mono{
   public:
	int xLoc;
	int yLoc;
	int wbc_roll;        //selectins
	int wbc_stick;       //integrens
	int wbc_migrate;     //diapedesis
	int mono_age;        //life span
	float TNFr;
	float IL_1r;
	float activation;
	int orientation;
	mono(int x, int y);
	mono(int x, int y, int age, int iTNFr, int iIL1r);
	void mono_function(int monoID);
	void mono_sniff();
};

struct mono_marrow{
   public:
	int xLoc;
	int yLoc;
	float wbc_roll;        //selectins
	float wbc_stick;       //integrens
	float wbc_migrate;     //diapedesis
	int pmn_age;         //life span
	float pmn_pcd;         //programmed cell death
	int mono_age;        //life span
	float TNFr;
	float IL_1r;
	float activation;
	mono_marrow(int x, int y);
	void mono_marrow_function();
	int orientation;
};

struct TH0{
   public:
	int xLoc;
	int yLoc;
	float wbc_roll;        //selectins
	float wbc_stick;       //integrens
	float wbc_migrate;     //diapedesis
	int pmn_age;         //life span
	float pmn_pcd;         //programmed cell death
	int mono_age;        //life span
	float TNFr;
	float IL_1r;
	int TH0_age;
	float activation;
	float proTH1;
	float proTH2;
	float rTH1;           //random holder for pro-TH1
	float rTH2;           //random holder for pro-TH2
	int orientation;
	TH0(int x, int y);
	TH0(int x, int y, int age);
	void TH0function(int index);
};

struct TH1{
   public:
	int xLoc;
	int yLoc;
	float wbc_roll;        //selectins
	float wbc_stick;       //integrens
	float wbc_migrate;     //diapedesis
	int pmn_age;         //life span
	float pmn_pcd;         //programmed cell death
	int mono_age;        //life span
	float TNFr;
	float IL_1r;
	int TH1_age;
	float activation;
	float pro_TH1;
	float pro_TH2;
	float rTH1;           //random holder for pro-TH1
	float rTH2;           //random holder for pro-TH2
	int orientation;
	TH1(int x, int y);
	TH1(int x, int y, int age);
	void TH1function(int index);
};

struct TH2{
   public:
	int xLoc;
	int yLoc;
	float wbc_roll;        //selectins
	float wbc_stick;       //integrens
	float wbc_migrate;     //diapedesis
	int pmn_age;         //life span
	float pmn_pcd;         //programmed cell death
	int mono_age;        //life span
	float TNFr;
	float IL_1r;
	int TH2_age;
	float activation;
	float pro_TH1;
	float pro_TH2;
	float rTH1;           //random holder for pro-TH1
	float rTH2;           //random holder for pro-TH2
	int orientation;
	TH2(int x, int y);
	TH2(int x, int y, int age);
	void TH2function(int index);
};

struct TH0_germ{
   public:
	int xLoc;
	int yLoc;
	int orientation;
	TH0_germ(int x, int y);
	void TH0_germ_function();
};

struct TH1_germ{
   public:
	int xLoc;
	int yLoc;
	int orientation;
	TH1_germ(int x, int y);
	void TH1_germ_function();
};

struct TH2_germ{
   public:
	int xLoc;
	int yLoc;
	int orientation;
	TH2_germ(int x, int y);
	void TH2_germ_function();
};

#endif
