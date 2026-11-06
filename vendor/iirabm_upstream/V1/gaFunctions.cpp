#include <vector>
#include "agents.h"
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include "Parameters.h"
#include "mpi.h"

using namespace std;


extern float PAFmult,TNFmult,sTNFrmult,IL1mult,sIL1rmult,IL1ramult,IFNgmult,IL4mult,IL6mult,
	IL8mult,IL10mult,IL12mult,GCSFmult;
extern const int tournamentSize,ga_popSize,procID;

extern mt19937 generator,gaGenerator;
extern uniform_int_distribution<int> distribution1000,distribution10,distribution9,distribution2,distributionBreed,distribution12,distributionMut;

extern vector<int> interventionLabels,breedIndex;
extern float possibleInterventions[];
extern float allInterventions[ga_popSize][genomeLength];
extern float allInterventionsReplicate[ga_popSize*replicateMult][genomeLength];
extern const float oxyCompare,infectCompare;
float breedPopulation[ga_popSize/tournamentSize][genomeLength];

float evalFitness(float oxy, float infect);
void fitnessSelection(float fitArray[]);
int nextGenCreation();
void crossover(float parent1[], float parent2[], float offspring[], int genomeLength);
int tournament(int contestants[tournamentSize],float fitArray[]);
void generateInterventions(float allInterventions[][genomeLength]);
void generateInterventionsRep(float allInterventions[][genomeLength],float allInterventionsReplicate[][genomeLength]);
float findNewValue(float in, int uord);
void mutate(float intervention[]);
int checkForConvergence();
void convertToReplicate(float allInterventions[][genomeLength],float allInterventionsReplicate[][genomeLength]);

float evalFitness(float oxy, float infect){
	float fitness;
	fitness=oxy+infect;
	return fitness;
}

void fitnessSelection(float fitArray[]){ //use "Tournament Selection" EOM 3.2.3

	int i,j,winner,x;
	int contestants[tournamentSize];	
	breedIndex.clear();
	shuffle(interventionLabels.begin(),interventionLabels.end(),generator);
	for(i=0;i<ga_popSize;i+=tournamentSize){	
//		cout<<"i="<<i<<"procID="<<procID<<"\n";
		for(j=0;j<tournamentSize;j++){
//			cout<<"j="<<j<<"procID="<<procID<<"\n";
			contestants[j]=interventionLabels[i+j];
		}
//		cout<<"Exit J Loop"<<"\n";
		winner=tournament(contestants,fitArray);
		breedIndex.push_back(winner);
	}
//	cout<<"Breed Size="<<breedIndex.size()<<"\n";
}

int nextGenCreation(){
	int i,j,size,r1,r2,answer=0,tempMut1,tempMut2;
	float parent1[genomeLength],parent2[genomeLength],tempOffspring1[genomeLength],
		tempOffspring2[genomeLength];
	
	answer=checkForConvergence();
	if(answer==1){
		return answer;
	}
	size=breedIndex.size();  //make sure to reset GA vectors
	shuffle(breedIndex.begin(),breedIndex.end(),generator);
	for(i=0;i<size;i++){
		for(j=0;j<genomeLength;j++){
			breedPopulation[i][j]=allInterventions[breedIndex[i]][j];		
		}
	}
	for(i=0;i<ga_popSize;i+=2){
		for(j=0;j<genomeLength;j++){
		r1=distributionBreed(generator);
		r2=distributionBreed(generator);
			parent1[j]=breedPopulation[r1][j]; 
			parent2[j]=breedPopulation[r2][j];
		}
		crossover(parent1,parent2,tempOffspring1,genomeLength);
		crossover(parent1,parent2,tempOffspring2,genomeLength);
		tempMut1=distribution10(generator);
		tempMut2=distribution10(generator);
		if(tempMut1<mutationRate){mutate(tempOffspring1);}
		if(tempMut2<mutationRate){mutate(tempOffspring2);}
		for(j=0;j<genomeLength;j++){
			allInterventions[i][j]=tempOffspring1[j];
			allInterventions[i+1][j]=tempOffspring2[j];
		}
	}
	return answer;
}

int tournament(int contestants[tournamentSize],float fitArray[ga_popSize]){
	float optimalFitness=999999;
	int i,answer;
	
	for(i=0;i<tournamentSize;i++){
		if(fitArray[contestants[i]]<optimalFitness){
			optimalFitness=fitArray[contestants[i]];
			answer=contestants[i];
		}
	}
	return answer;
}

void crossover(float parent1[], float parent2[], float offspring[], int genomeLength){
	int i,temp;
	
	for(i=0;i<genomeLength;i++){
		temp=distribution2(generator);
		if(temp==0){
			offspring[i]=parent1[i];
		} else {
			offspring[i]=parent2[i];
		}
	}
}

void mutate(float intervention[]){
	int i,tempIndex,upordown,mutValue;
	float mutation;
	tempIndex=distributionMut(generator);
	upordown=distribution2(generator);
      mutation=findNewValue(intervention[tempIndex],upordown);
      intervention[tempIndex]=mutation;
//	intervention[tempIndex]=possibleInterventions[mutValue];
}

void generateInterventions(float allInterventions[][genomeLength]){
	int i,j,temp;
	for(i=0;i<ga_popSize;i++){
		for(j=0;j<genomeLength;j++){
			temp=distribution9(gaGenerator);
			allInterventions[i][j]=possibleInterventions[temp];	
		}
	}
}

void generateInterventionsRep(float allInterventions[][genomeLength],float allInterventionsReplicate[][genomeLength]){
	int i,j,k,temp;
	for(i=0;i<ga_popSize;i++){
		for(j=0;j<genomeLength;j++){
			temp=distribution9(gaGenerator);
			allInterventions[i][j]=possibleInterventions[temp];	
			for(k=0;k<replicateMult;k++){
				allInterventionsReplicate[i*replicateMult+k][j]=possibleInterventions[temp];
			}	
		}
	}
}

float findNewValue(float in, int uord){
	int i,position;
	for(i=0;i<=8;i++){
		if(in==possibleInterventions[i]){
			position=i;
			break;
		}
	}
	if(uord==0){position--;}
	if(uord==1){position++;}
	if(position==-1){position=0;}
	if(position==9){position=8;}
	return possibleInterventions[position];
}

int checkForConvergence(){
	int i,j,answer=1;
	float temp;
	
	for(i=0;i<genomeLength;i++){
		temp=allInterventions[0][i];
		for(j=1;j<breedIndex.size();j++){
			if(allInterventions[breedIndex[j]][i]!=temp){
				answer=0;
				return answer;
			}
		}
	}
	return answer;
}

void convertToReplicate(float allInterventions[][genomeLength],float allInterventionsReplicate[][genomeLength]){
	int i,i2,j;
	float tempInt[genomeLength];
	for(i=0;i<ga_popSize;i++){
		for(j=0;j<genomeLength;j++){
			tempInt[j]=allInterventions[i][j];
		}
		for(i2=i*replicateMult;i2<i+replicateMult;i2++){
			for(j=0;j<genomeLength;j++){
				allInterventionsReplicate[i2][j]=tempInt[j];
			}
		}
	}

}

