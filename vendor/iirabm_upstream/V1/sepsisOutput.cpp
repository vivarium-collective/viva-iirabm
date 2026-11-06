#include <vector>
#include "agents.h"
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include "Parameters.h"

extern ofstream ecMap,oxyInfectFinalFile,dailyOutputFile,allOutputFile,gaFile;
extern vector<EC> ecArray;
extern vector<int> breedIndex;
extern int cellGrid[101][101],procID;
extern float allInterventions[ga_popSize][genomeLength];
extern float fitnessValues[ga_popSize];
void infectionMap(int time);
void outputSetup();
void outputFinalize();

void oxy_infectFinalWriter(int seed, int inj_number, int infect_spread, int numRecurInj, 
	int numInfectRepeat, float oxyheal, float oxy, float infect){
	oxyInfectFinalFile<<seed<<" "<<inj_number<<" "<<infect_spread<<" "<<numRecurInj<<" "
	<<numInfectRepeat<<" "<<oxyheal<<" "<<oxy<<" "<<infect<<"\n";
}

void infectionMap(int time){
	int i,size;
		
	size=ecArray.size();
	cout<<"ProcID="<<procID<<" size="<<size<<"\n";
	for(i=0;i<size;i++){
		ecMap<<time<<" "<<ecArray[i].xLoc<<" "<<ecArray[i].yLoc<<" "<<ecArray[i].infection<<
			" "<<ecArray[i].oxy<<" "<<ecArray[i].ec_activation<<" "<<ecArray[i].ec_roll
			<<" "<<ecArray[i].ec_stick<<" "<<ecArray[i].ec_migrate<<" "<<ecArray[i].cytotox
			<<" "<<ecArray[i].endotoxin<<" "<<ecArray[i].PAF<<" "<<ecArray[i].TNF
			<<" "<<ecArray[i].sTNFr<<" "<<ecArray[i].IL1<<" "<<ecArray[i].sIL1r
			<<" "<<ecArray[i].IL1ra<<" "<<ecArray[i].IFNg
			<<" "<<ecArray[i].IL4<<" "<<ecArray[i].IL8
			<<" "<<ecArray[i].IL10<<" "<<ecArray[i].IL12<<" "<<ecArray[i].GCSF<<"\n";;
	}	
}

void dailyOutputWriter(int* intRecv,float* flRecv,int* pmnPop,int* monoPop,int* TH1pop,int* TH2pop,
	float* dTNF,float* dsTNFr,float* dIL10,float* dGCSF,float* dpTH1,float* dpTH2,
	float* dOxy, float* dInfect, float* dPAF, float* dIFNg, float* dIL4, float* dIL8,
	float* dIL12, float* dsIL1r, float* dIL1ra, float* dIL1, int numProcs){
	
	int i,j;
	
	for(i=0;i<numProcs;i++){
		dailyOutputFile<<intRecv[i*5]<<" "<<intRecv[i*5+1]<<" "<<intRecv[i*5+2]<<" "<<intRecv[i*5+3]<<
			" "<<intRecv[i*5+4]<<" "<<flRecv[i]<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<pmnPop[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<monoPop[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<TH1pop[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<TH2pop[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		///////////////////////////
		for(j=0;j<112;j++){
			dailyOutputFile<<dOxy[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dInfect[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dTNF[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dsTNFr[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dIL10[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dGCSF[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dpTH1[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dpTH2[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dPAF[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dIFNg[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dIL4[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dIL8[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dIL12[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dsIL1r[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dIL1ra[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
		for(j=0;j<112;j++){
			dailyOutputFile<<dIL1[i*112+j]<<" ";
		}
		dailyOutputFile<<"\n";
	}
}

void allOutputWriter(int* intRecv,float* flRecv,int* pmnPop,int* monoPop,int* TH1pop,int* TH2pop,
	float* dTNF,float* dsTNFr,float* dIL10,float* dGCSF,float* dpTH1,float* dpTH2,
	float* dOxy, float* dInfect, float* dPAF, float* dIFNg, float* dIL4, float* dIL8,
	float* dIL12, float* dsIL1r, float* dIL1ra, float* dIL1, int numProcs){
	
	int numTimeSteps=numWriteSteps;
	
	int i,j;
	
	for(i=0;i<numProcs;i++){
		allOutputFile<<intRecv[i*5]<<" "<<intRecv[i*5+1]<<" "<<intRecv[i*5+2]<<" "<<intRecv[i*5+3]<<
			" "<<intRecv[i*5+4]<<" "<<flRecv[i]<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<pmnPop[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<monoPop[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<TH1pop[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<TH2pop[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		///////////////////////////
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dOxy[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dInfect[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dTNF[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dsTNFr[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dIL10[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dGCSF[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dpTH1[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dpTH2[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dPAF[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dIFNg[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dIL4[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dIL8[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dIL12[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dsIL1r[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dIL1ra[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
		for(j=0;j<numTimeSteps;j++){
			allOutputFile<<dIL1[i*numTimeSteps+j]<<" ";
		}
		allOutputFile<<"\n";
	}
}

void GAwriter(float fitArray[]){
	int i,j;
	for(i=0;i<breedIndex.size();i++){
		gaFile<<fitArray[breedIndex[i]]<<" ";
		for(j=0;j<genomeLength;j++){
			gaFile<<allInterventions[breedIndex[i]][j]<<" ";
		}
		gaFile<<"\n";
	}
}

void movieSetup(){
	char fileExt[]=".dat";
	char fileStart[]="modelState_";
	char buffer[256];
	
	sprintf(buffer, "%04d", procID);
	string str(buffer);
	strcat(fileStart,buffer);
	strcat(fileStart,fileExt);
	string fileName(fileStart);
	ecMap.open(fileName.c_str());
}

void outputSetup(){
	char fileExt[]=".dat";
	char fileStart[]="modelState2_";
	char buffer[256];
	
	sprintf(buffer, "%04d", procID);
	string str(buffer);
	strcat(fileStart,buffer);
	strcat(fileStart,fileExt);
	if(oxyInfectFinal==1){
		oxyInfectFinalFile.open(fileName);}
	if(dailyOutput==1){
		dailyOutputFile.open(fileName);}
	if(writeEverything==1){
			allOutputFile.open(fileName);}
	string fileName(fileStart);
	if(movie==1){
		ecMap.open(fileName.c_str());}
	if(GA==1){
		gaFile.open(fileNameGA.c_str());}
}

void outputFinalize(){
	if(oxyInfectFinal==1){
		oxyInfectFinalFile.close();}
	if(writeEverything==1){
		allOutputFile.close();}
	if(movie==1){
		ecMap.close();}
}

void outputClose(){
	if(oxyInfectFinal==1){
		oxyInfectFinalFile.close();}
	if(dailyOutput==1){
		dailyOutputFile.close();}
	if(writeEverything==1){
			allOutputFile.close();}
	if(movie==1){
		ecMap.close();}
	if(GA==1){
		gaFile.close();}
}
