#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include "mpi.h"
#include "agents.h"
#include "Parameters.h"
using namespace std;

extern const int cellCapacity,xDim,yDim,injuryStep,parameterInput,numTimeSteps,numWriteSteps,
	ga_popSize,maxIter;
extern int numProcs;
extern const float antibioticMultiplier;

vector<EC> ecArray; 
vector<int> ecIndexes;     
vector<pmn> pmnArray;
vector<mono> monoArray;
vector<TH0> TH0array;
vector<TH1> TH1array;
vector<TH2> TH2array;
vector<pmn_marrow> pmn_marrowArray;
vector<mono_marrow> mono_marrowArray;
vector<TH0_germ> TH0_germArray;
vector<TH1_germ> TH1_germArray;
vector<TH2_germ> TH2_germArray;
vector<int> interventionLabels,breedIndex;
int procID,cellGrid[101][101],p1[10000],p3[10000],p4[10000],p5[10000],p6[10000],np1,numProcs;
float p2[10000];
int* np1ptr;

float system_oxy,oxyDeficit,totalInfection,total_TNF,total_sTNFr,total_IL10,
      total_IL6,total_GCSF,total_proTH1,total_proTH2,total_IFNg,total_PAF,
      total_IL1,total_IL4,total_IL8,total_IL12,total_sIL1r,total_IL1ra;
float PAFmult,TNFmult,sTNFrmult,IL1mult,sIL1rmult,IL1ramult,IFNgmult,IL4mult,
	IL8mult,IL10mult,IL12mult,GCSFmult;
float possibleInterventions[]={0,0.25,0.5,0.75,1,2,5,10,20};	
float d_oxy[numWriteSteps],d_infect[numWriteSteps],d_TNF[numWriteSteps],d_sTNFr[numWriteSteps],
	d_IL10[numWriteSteps],d_IL6[numWriteSteps],d_GCSF[numWriteSteps],d_pTH1[numWriteSteps],
	d_pTH2[numWriteSteps],d_IFNg[numWriteSteps],d_PAF[numWriteSteps],d_IL1[numWriteSteps],
	d_IL4[numWriteSteps],d_IL8[numWriteSteps],d_IL12[numWriteSteps],d_IL1ra[numWriteSteps],
	d_sIL1r[numWriteSteps];
int d_pmn[numWriteSteps],d_mono[numWriteSteps],d_TH1[numWriteSteps],d_TH2[numWriteSteps];
float allInterventions[ga_popSize][genomeLength];
float allInterventionsReplicate[ga_popSize*replicateMult][genomeLength];
float myIntervention[genomeLength];
float testInt[testLength];

float mGCSF,mPAF,mTNF,mSTNFR,mIL1,mSIL1R,mIL1RA,mIFNg,mIL4,mIL8,mIL10,mIL12;

int reseedNum;

ofstream ecMap,oxyInfectFinalFile,dailyOutputFile,allOutputFile,gaFile;
mt19937 generator,gaGenerator;
uniform_int_distribution<int> distribution10k(0,9999);
uniform_int_distribution<int> distribution1000(0,999);
uniform_int_distribution<int> distribution100(0,99);
uniform_int_distribution<int> distribution50(0,49);
uniform_int_distribution<int> distribution12(0,11);
uniform_int_distribution<int> distribution10(0,9);
uniform_int_distribution<int> distribution9(0,8);
uniform_int_distribution<int> distribution8(0,7);
uniform_int_distribution<int> distribution5(0,4);
uniform_int_distribution<int> distribution3(0,2);
uniform_int_distribution<int> distribution2(0,1);
uniform_int_distribution<int> distributionMut(0,genomeLength-1);
uniform_int_distribution<int> distributionBreed(0,ga_popSize/tournamentSize-1);

void mainSimulation(int inj_number, float oxyHeal, int infectSpread, int numRecurInj, int numInfectRepeat, int seed);
void postSimulationData(int index);
void printIntervention();

extern void updateSystemOxy(int step);
extern void evaporate();
extern void diffuse();
extern void initialize(); 
extern void injure_infectionFRD(int inj_number);
extern void infectionMap(int time);
extern void oxy_infectFinalWriter(int seed, int inj_number, int infect_spread, int numRecurInj, 
	int numInfectRepeat, float oxyheal, float oxy, float infect);
extern void pickParams(int p1[], int* np1, float p2[], int p3[], int p4[], int p5[], int p6[], 
	int procID, int maxProcs);
extern void recur_injury();
extern void dailyOutputWriter(int* intRecv,float* flRecv,int* pmnPop,int* monoPop,int* TH1pop,int* TH2pop,
	float* dTNF,float* dsTNFr,float* dIL10,float* dGCSF,float* dpTH1,float* dpTH2,
	float* dOxy, float* dInfect, float* dPAF, float* dIFNg, float* dIL4, float* dIL8,
	float* dIL12, float* dsIL1ra, float* dIL1ra, float* dIL1, int numProcs);
extern void allOutputWriter(int* intRecv,float* flRecv,int* pmnPop,int* monoPop,int* TH1pop,int* TH2pop,
	float* dTNF,float* dsTNFr,float* dIL10,float* dGCSF,float* dpTH1,float* dpTH2,
	float* dOxy, float* dInfect, float* dPAF, float* dIFNg, float* dIL4, float* dIL8,
	float* dIL12, float* dsIL1ra, float* dIL1ra, float* dIL1, int numProcs);
extern void readParams(int procID, int nProc,int p1[], int* np1, float p2[], int p3[], int p4[], 
	int p5[], int p6[]);
extern void applyAntibiotics();
extern void outputSetup();
extern void outputFinalize();
extern void generateInterventions(float allInterventions[][genomeLength]);
extern void generateInterventionsRep(float allInterventions[][genomeLength],float allInterventionsReplicate[][genomeLength]);
extern float evalFitness(float oxy, float infect);
extern int nextGenCreation();
extern void fitnessSelection(float fitArray[]);
extern void GAwriter(float fitArray[]);
extern void spreadParams(int p1[], float p2[], int p3[], int p4[], int p5[], int p6[], int* np1);
extern void outputClose();
extern int checkForConvergence();
extern void applyIntervention(int flag);
extern void applyTestIntervention(int flag);
extern void movieSetup();
extern void clearIntervention();
extern void convertToReplicate(float allInterventions[][genomeLength],float allInterventionsReplicate[][genomeLength]);
extern void readIntervention();

int main(){ 	
	int i1,i,iI,numGA_Iter=0,ifit,jfit,kfit;
	int gaFlag=0,x,rseed;
	float myFitness=0;
	double starttime,endtime;
	float *fitArray;
	np1=0;
	np1ptr=&np1;

	MPI::Init();
  	procID=MPI::COMM_WORLD.Get_rank();
  	numProcs=MPI::COMM_WORLD.Get_size();
  	
  	starttime=MPI::Wtime();
  	
  	if(procID==0){
  		if(movie==0){outputSetup();}
  		cout<<"numProcs="<<numProcs<<"\n";
		for(i=0;i<ga_popSize;i++){
			interventionLabels.push_back(i);
		}
  	}
  	if(movie==1){
  		movieSetup();
  	}
  	
  	if(testIntervention==1){
  		readIntervention();
  	}
  	
  	if(GA==0&&numReseeds>0){
  		spreadParams(p1,p2,p3,p4,p5,p6,np1ptr);
  		mainSimulation(p1[0],p2[0],p3[0],p4[0],p5[0],p6[0]);
  		postSimulationData(0);
  	}

	if(GA==0&&numReseeds==0){
		if(parameterInput==0){
			pickParams(p1,np1ptr,p2,p3,p4,p5,p6,procID,numProcs);}
		else{
			readParams(procID,numProcs,p1,np1ptr,p2,p3,p4,p5,p6);}
		for(i1=0;i1<np1;i1++){
			cout<<"ProcID="<<procID<<" Params="<<p1[i1]<<" "<<p2[i1]<<" "<<p3[i1]<<" "<<p4[i1]<<" "<<p5[i1]<<" "<<p6[i1]<<"\n";
			mainSimulation(p1[i1],p2[i1],p3[i1],p4[i1],p5[i1],p6[i1]);
			cout<<"oxy="<<oxyDeficit<<"\n";
                  cout<<"Infect="<<totalInfection<<"\n";
			postSimulationData(i1);
		}
	}
	if(GA==1){
		gaGenerator.seed(GAseed);
		spreadParams(p1,p2,p3,p4,p5,p6,np1ptr);
		if(procID==0){
			generateInterventions(allInterventions);	
			fitArray=new float[ga_popSize];	
		}

		MPI::COMM_WORLD.Scatter(allInterventions,genomeLength,MPI::FLOAT,myIntervention,genomeLength,MPI::FLOAT,0);

		while((gaFlag==0)&&(numGA_Iter<maxIter)){
			myFitness=0;
			if(reseedPostInt==0){
				for(rseed=seedStart;rseed<=seedMax;rseed+=seedIncrement){
					mainSimulation(p1[0],p2[0],p3[0],p4[0],p5[0],rseed);
					if(dailyOutput==1||writeEverything==1){
						postSimulationData(0);
					}
					myFitness+=evalFitness(oxyDeficit,totalInfection);
				}
			}
			if(reseedPostInt>0){
				for(reseedNum=1;reseedNum<=replicateMult;reseedNum++){
					rseed=seedStart;
					mainSimulation(p1[0],p2[0],p3[0],p4[0],p5[0],rseed);
					myFitness+=evalFitness(oxyDeficit,totalInfection);
					cout<<"ProcID="<<procID<<" fit="<<myFitness<<" oxy="<<oxyDeficit<<" TI="<<totalInfection<<"\n";
				}
			}
			MPI::COMM_WORLD.Gather(&myFitness,1,MPI::FLOAT,fitArray,1,MPI::FLOAT,0);

			if(procID==0){
				fitnessSelection(fitArray);
				GAwriter(fitArray);
				gaFlag=nextGenCreation();
			}
			
			MPI::COMM_WORLD.Scatter(allInterventions,genomeLength,MPI::FLOAT,myIntervention,genomeLength,MPI::FLOAT,0);
		
			MPI::COMM_WORLD.Bcast(&gaFlag,1,MPI::INT,0);
			numGA_Iter++;
			if(procID==0){cout<<"Number of Iterations=  "<<numGA_Iter<<"\n";}
		}
		if(procID==0){
			delete [] fitArray;
		}	
	}

	if(procID==0){
		endtime=MPI::Wtime();
		cout<<"Time="<<endtime-starttime<<"\n";
 		cout<<"All Simulations Complete\n";	
 		outputClose();
 	}
 	MPI::Finalize();	
};

void mainSimulation(int inj_number, float oxyHeal, int infectSpread, int numRecurInj, int numInfectRepeat, int seed){

      int length,i,j,step,ii,iend,dOut,count6hr,antibiotic1,antibiotic2,numAntibiotics,aFlag,iFlag,istep;
      int iFlag2,mstep2,iFlagTest;
      float fitness;
	generator.seed(seed);
	cout<<distribution100(generator)<<"\n";
	initialize();
	dOut=0;
	step=0;
	istep=0;
	count6hr=0;
	antibiotic1=0;
	antibiotic2=0;
	numAntibiotics=0;
	aFlag=0;
	iFlag=0;
	iFlag2=0;
	iFlagTest=0;
	mstep2=movieStep;
	for(i=0;i<xDim*yDim;i++){
		ecIndexes.push_back(i);
	}
	for(i=0;i<numTimeSteps;i++){
		if((numReseeds>0)&&(istep==reseedTime)&&(procID!=0)){
			generator.seed(procID+1);
		}
		if((reseedPostInt>0)&&(istep==reseedTime)&&(GA==1)){
			generator.seed(reseedNum);
		}

		if(step==0){injure_infectionFRD(inj_number);}
		step++;
		istep++;
		if(step==injuryStep){step=1;}
		antibiotic1++;
		antibiotic2++;

		updateSystemOxy(istep);
		if(oxyDeficit>8161){
			if((writeEverything==1)||(dailyOutput==1)){
			for(iend=dOut;iend<numWriteSteps;iend++){
				d_oxy[iend]=-1.0;
				d_infect[iend]=-1.0;
				d_TNF[iend]=-1.0;
				d_sTNFr[iend]=-1.0;
				d_IL10[iend]=-1.0;
				d_GCSF[iend]=-1.0;
				d_pTH1[iend]=-1.0;
				d_pTH2[iend]=-1.0;
				d_pmn[iend]=-1;
				d_mono[iend]=-1;
				d_TH1[iend]=-1;
				d_TH2[iend]=-1;
				d_IFNg[iend]=-1.0;
				d_PAF[iend]=-1.0;
				d_IL1[iend]=-1.0;
				d_IL4[iend]=-1.0;
				d_IL8[iend]=-1.0;
				d_IL12[iend]=-1.0;
				d_sIL1r[iend]=-1.0;
				d_IL1ra[iend]=-1.0;

			}
			}
			return;
		}	

		if((movie==1)&&(istep==mstep2)){
			infectionMap(istep);
			mstep2=mstep2+51;
		}

		if(((dailyOutput==1)&&(count6hr==50))||(writeEverything==1)){
				d_oxy[dOut]=oxyDeficit;
				d_infect[dOut]=totalInfection;
				d_TNF[dOut]=total_TNF;
				d_sTNFr[dOut]=total_sTNFr;
				d_IL10[dOut]=total_IL10;
				d_GCSF[dOut]=total_GCSF;
				d_pTH1[dOut]=total_proTH1;
				d_pTH2[dOut]=total_proTH2;
				d_pmn[dOut]=pmnArray.size();
				d_mono[dOut]=monoArray.size();
				d_TH1[dOut]=TH1array.size();
				d_TH2[dOut]=TH2array.size();
				d_IFNg[dOut]=total_IFNg;
				d_PAF[dOut]=total_PAF;
				d_IL1[dOut]=total_IL1;
				d_IL4[dOut]=total_IL4;
				d_IL8[dOut]=total_IL8;
				d_IL12[dOut]=total_IL12;
				d_sIL1r[dOut]=total_sIL1r;
				d_IL1ra[dOut]=total_IL1ra;
				dOut++;
		}

		length=TH0array.size();
 		if(length>0){
 			shuffle(TH0array.begin(),TH0array.end(),generator);}
		j=0;
		while(j<length){
			TH0array[j].TH0function(j);
			j++;
			length=TH0array.size();}
		length=ecArray.size();
		shuffle(ecIndexes.begin(),ecIndexes.end(),generator);
		j=0;
		while(j<length){
			ecArray[ecIndexes[j]].inj_function(infectSpread,numInfectRepeat);
			ecArray[ecIndexes[j]].ECfunction(oxyHeal);
			j++;
			length=ecArray.size();}	
		length=pmnArray.size();
 		if(length>0){
 			shuffle(pmnArray.begin(),pmnArray.end(),generator);}
		j=0;
		while(j<length){
			pmnArray[j].pmn_function(j);
			j++;
			length=pmnArray.size();}
		length=monoArray.size();
  		if(length>0){
  			shuffle(monoArray.begin(),monoArray.end(),generator);}
		j=0;
		while(j<length){
			monoArray[j].mono_function(j);
			j++;
			length=monoArray.size();}			
		length=TH1array.size();
  		if(length>0){
  			shuffle(TH1array.begin(),TH1array.end(),generator);}
		j=0;
		while(j<length){
			TH1array[j].TH1function(j);
			j++;
			length=TH1array.size();}
		length=TH2array.size();
  		if(length>0){
  			shuffle(TH2array.begin(),TH2array.end(),generator);}
		j=0;
		while(j<length){
			TH2array[j].TH2function(j);
			j++;
			length=TH2array.size();}
		length=pmn_marrowArray.size();
  		if(length>0){
  			shuffle(pmn_marrowArray.begin(),pmn_marrowArray.end(),generator);}
		j=0;
		while(j<length){
			pmn_marrowArray[j].pmn_marrow_function();
			j++;
			length=pmn_marrowArray.size();}	
		
		length=mono_marrowArray.size();
  		if(length>0){
  			shuffle(mono_marrowArray.begin(),mono_marrowArray.end(),generator);}
		j=0;
		while(j<length){
			mono_marrowArray[j].mono_marrow_function();
			j++;
			length=mono_marrowArray.size();}
			
		length=TH1_germArray.size();
  		if(length>0){
  			shuffle(TH1_germArray.begin(),TH1_germArray.end(),generator);}
		j=0;
		while(j<length){
			TH1_germArray[j].TH1_germ_function();
			j++;
			length=TH1_germArray.size();}
			
		length=TH2_germArray.size();
  		if(length>0){
  			shuffle(TH2_germArray.begin(),TH2_germArray.end(),generator);}
		j=0;
		while(j<length){
			TH2_germArray[j].TH2_germ_function();
			j++;
			length=TH2_germArray.size();}
		
		length=TH0_germArray.size();
  		if(length>0){
  			shuffle(TH0_germArray.begin(),TH0_germArray.end(),generator);}
		j=0;
		while(j<length){
			TH0_germArray[j].TH0_germ_function();
			j++;
			length=TH0_germArray.size();}
			
			
		evaporate();		

		if(step==injuryStep-1){
			for(ii=1;ii<=numRecurInj;ii++){
				recur_injury();
			}
		}
		
		diffuse();
		
		if(antibioticMultiplier>0){
			if((antibiotic1==102)&&(aFlag==0)){
				aFlag=1;
				applyAntibiotics();
				antibiotic2=0;
				}
			if((antibiotic2==injuryStep)&&(aFlag==1)&&(numAntibiotics<1100)){
				numAntibiotics++;
				applyAntibiotics();
				antibiotic2=0;
			}
		}
		
		if(GA==1&&testIntervention==0){
			if((istep==interventionTimeStart+iFlag*interventionTimeSeperation)&&(iFlag<numIndependentInterventions)){
				iFlag++;
				applyIntervention(iFlag-1);
//				printIntervention();				
			}
			if((istep==interventionTimeStart+interventionLength+(iFlag-1)*interventionTimeSeperation)){
				clearIntervention();
			}
		}

		if(testIntervention==1&&GA==0){
			if((istep==interventionTimeStart+iFlag*interventionTimeSeperation)&&(iFlag<numIndependentInterventions)){
				applyTestIntervention(iFlag);
				printIntervention();
				iFlag++;
			}
			if((istep==interventionTimeStart+interventionLength+(iFlag-1)*interventionTimeSeperation)){
				clearIntervention();
			}
		}
		if(testIntervention==1&&GA==1){
			if((istep==testStart+iFlag*interventionTimeSeperation)&&(iFlagTest<numTests)){
				applyTestIntervention(iFlag);
				iFlagTest++;
			}
			if((istep==testStart+interventionLength+(iFlagTest-1)*interventionTimeSeperation)){
				clearIntervention();
			}
			
			if((istep==interventionTimeStart+iFlag*interventionTimeSeperation)&&(iFlag<numIndependentInterventions)){
				applyIntervention(iFlag);
				printIntervention();
				iFlag++;
			}
			if((istep==interventionTimeStart+interventionLength+(iFlag-1)*interventionTimeSeperation)){
				clearIntervention();
			}
			
		}
		
	}
}

void postSimulationData(int index){
	int intSend[5],i;
	float flSend[3],healSend[1];
	int *intRecv,*pmnPop,*monoPop,*TH1pop,*TH2pop;
	float *flRecv,*dOxy,*dInfect,*dTNF,*dsTNFr,*dIL10,*dIL6,*dGCSF,*dPTH1,*dPTH2,*dPAF,*dIFNg,
		*dIL1,*dIL4,*dIL8,*dIL12,*dsIL1r,*dIL1ra;

	if(oxyInfectFinal==1){
		if(procID==0){
			intRecv=new int[5*numProcs];
			flRecv=new float[3*numProcs];
		}
	}

	if((dailyOutput==1)||(writeEverything==1)){
		if(procID==0){
			intRecv=new int[5*numProcs];
			flRecv=new float[numProcs];
			pmnPop=new int[numWriteSteps*numProcs];
			monoPop=new int[numWriteSteps*numProcs];
			TH1pop=new int[numWriteSteps*numProcs];
			TH2pop=new int[numWriteSteps*numProcs];
			dOxy=new float[numWriteSteps*numProcs];
			dInfect=new float[numWriteSteps*numProcs];
			dTNF=new float[numWriteSteps*numProcs];
			dsTNFr=new float[numWriteSteps*numProcs];
			dIL10=new float[numWriteSteps*numProcs];
			dIL6=new float[numWriteSteps*numProcs];
			dGCSF=new float[numWriteSteps*numProcs];
			dPTH1=new float[numWriteSteps*numProcs];
			dPTH2=new float[numWriteSteps*numProcs];
			dPAF=new float[numWriteSteps*numProcs];
			dIFNg=new float[numWriteSteps*numProcs];
			dsIL1r=new float[numWriteSteps*numProcs];
			dIL1ra=new float[numWriteSteps*numProcs];
			dIL1=new float[numWriteSteps*numProcs];
			dIL4=new float[numWriteSteps*numProcs];
			dIL8=new float[numWriteSteps*numProcs];
			dIL12=new float[numWriteSteps*numProcs];
		}
	}

	intSend[0]=p6[index];
	intSend[1]=p1[index];
	intSend[2]=p3[index];
	intSend[3]=p4[index];
	intSend[4]=p5[index];
 		
	if(oxyInfectFinal==1){
 		flSend[0]=p2[index];
 		flSend[1]=oxyDeficit;
 		flSend[2]=totalInfection;
		MPI::COMM_WORLD.Gather(intSend,5,MPI::INT,intRecv,5,MPI::INT,0);
		MPI::COMM_WORLD.Gather(flSend,3,MPI::FLOAT,flRecv,3,MPI::FLOAT,0);
		if(procID==0){
			for(i=0;i<numProcs;i++){
				oxy_infectFinalWriter(intRecv[i*5],intRecv[i*5+1],intRecv[i*5+2],intRecv[i*5+3],
					intRecv[i*5+4],flRecv[i*3],flRecv[i*3+1],flRecv[i*3+2]);							
		}}
		oxyInfectFinalFile.flush();
	}
	if((dailyOutput==1)||(writeEverything==1)){
		healSend[0]=p2[index];
		MPI::COMM_WORLD.Gather(intSend,5,MPI::INT,intRecv,5,MPI::INT,0);
		MPI::COMM_WORLD.Gather(healSend,1,MPI::FLOAT,flRecv,1,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_pmn,numWriteSteps,MPI::INT,pmnPop,numWriteSteps,MPI::INT,0);
		MPI::COMM_WORLD.Gather(d_mono,numWriteSteps,MPI::INT,monoPop,numWriteSteps,MPI::INT,0);
		MPI::COMM_WORLD.Gather(d_TH1,numWriteSteps,MPI::INT,TH1pop,numWriteSteps,MPI::INT,0);
		MPI::COMM_WORLD.Gather(d_TH2,numWriteSteps,MPI::INT,TH2pop,numWriteSteps,MPI::INT,0);
		MPI::COMM_WORLD.Gather(d_TNF,numWriteSteps,MPI::FLOAT,dTNF,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_sTNFr,numWriteSteps,MPI::FLOAT,dsTNFr,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IL10,numWriteSteps,MPI::FLOAT,dIL10,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IL6,numWriteSteps,MPI::FLOAT,dIL6,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_GCSF,numWriteSteps,MPI::FLOAT,dGCSF,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_pTH1,numWriteSteps,MPI::FLOAT,dPTH1,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_pTH2,numWriteSteps,MPI::FLOAT,dPTH2,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_oxy,numWriteSteps,MPI::FLOAT,dOxy,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_infect,numWriteSteps,MPI::FLOAT,dInfect,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_PAF,numWriteSteps,MPI::FLOAT,dPAF,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IFNg,numWriteSteps,MPI::FLOAT,dIFNg,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_sIL1r,numWriteSteps,MPI::FLOAT,dsIL1r,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IL1ra,numWriteSteps,MPI::FLOAT,dIL1ra,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IL1,numWriteSteps,MPI::FLOAT,dIL1,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IL4,numWriteSteps,MPI::FLOAT,dIL4,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IL8,numWriteSteps,MPI::FLOAT,dIL8,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Gather(d_IL12,numWriteSteps,MPI::FLOAT,dIL12,numWriteSteps,MPI::FLOAT,0);
		MPI::COMM_WORLD.Barrier();
		if(procID==0){
			if(dailyOutput==1){
				dailyOutputWriter(intRecv,flRecv,pmnPop,monoPop,TH1pop,TH2pop,
					dTNF,dsTNFr,dIL10,dGCSF,dPTH1,dPTH2,dOxy,dInfect,
					dPAF,dIFNg,dIL4,dIL8,dIL12,dsIL1r,dIL1ra,dIL1,numProcs);
			}
			if(writeEverything==1){
				allOutputWriter(intRecv,flRecv,pmnPop,monoPop,TH1pop,TH2pop,
					dTNF,dsTNFr,dIL10,dGCSF,dPTH1,dPTH2,dOxy,dInfect,
					dPAF,dIFNg,dIL4,dIL8,dIL12,dsIL1r,dIL1ra,dIL1,numProcs);
	}}}
	if((dailyOutput==1)||(writeEverything==1)){	
		if(procID==0){
			delete [] dOxy;
			delete [] dInfect;
			delete [] dTNF;
			delete [] dsTNFr;
			delete [] dIL10;
			delete [] dIL6;
			delete [] dGCSF;
			delete [] dPTH1;
			delete [] dPTH2;
			delete [] dPAF;
			delete [] dIFNg;
			delete [] dsIL1r;
			delete [] dIL1ra;
			delete [] dIL1;
			delete [] dIL4;
			delete [] dIL8;
			delete [] dIL12;
			delete [] intRecv;
			delete [] flRecv;
			delete [] pmnPop;
			delete [] monoPop;
			delete [] TH1pop;
			delete [] TH2pop;
		}
	}
}

void printIntervention(){

	cout<<"ProcID="<<procID<<" PAF="<<PAFmult<<"\n";
	cout<<"ProcID="<<procID<<" TNF="<<TNFmult<<"\n";
	cout<<"ProcID="<<procID<<" STNFR="<<sTNFrmult<<"\n";
	cout<<"ProcID="<<procID<<" IL1="<<IL1mult<<"\n";
	cout<<"ProcID="<<procID<<" SIL1R="<<sIL1rmult<<"\n";
	cout<<"ProcID="<<procID<<" IL1RA="<<IL1ramult<<"\n";
	cout<<"ProcID="<<procID<<" IFNg="<<IFNgmult<<"\n";
	cout<<"ProcID="<<procID<<" IL4="<<IL4mult<<"\n";
	cout<<"ProcID="<<procID<<" IL8="<<IL8mult<<"\n";
	cout<<"ProcID="<<procID<<" Il10="<<IL10mult<<"\n";
	cout<<"ProcID="<<procID<<" IL12="<<IL12mult<<"\n";
	cout<<"ProcID="<<procID<<" GCSF="<<GCSFmult<<"\n";	
}

