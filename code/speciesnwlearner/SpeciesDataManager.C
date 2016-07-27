
#include <string.h>
#include "Variable.H"
#include "Evidence.H"
#include "Error.H"
#include "VariableManager.H"
#include "Potential.H"
#include "PotentialManager.H"
#include "SlimFactor.H"
#include "FactorGraph.H"
#include "EvidenceManager.H"
#include "SpeciesDataManager.H"


SpeciesDataManager::SpeciesDataManager()
{
}

SpeciesDataManager::~SpeciesDataManager()
{
}

int
SpeciesDataManager::setVariableManager(VariableManager* aPtr)
{
	varMgr=aPtr;
	return 0;
}

int
SpeciesDataManager::createFactorGraph()
{
	VSET& variableSet=varMgr->getVariableSet();
	int globalFactorID=0;
	fgraph=new FactorGraph;
	//Create the factors and the Markov blanket variables using the neighbours of each variable
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		SlimFactor* sFactor=new SlimFactor;
		sFactor->vIds=new int[1];
		sFactor->vIds[0]=vIter->first;
		sFactor->vCnt=1;
		sFactor->fId=globalFactorID;
		globalFactorID++;
		fgraph->setFactor(sFactor);
	}
	return 0;
}

int
SpeciesDataManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}


int 
SpeciesDataManager::setPotentialManager(PotentialManager* aPtr)
{
	potMgr=aPtr;
	return 0;
}

int
SpeciesDataManager::setOutputLoc(const char* aPtr)
{
	strcpy(outputLoc,aPtr);
	return 0;
}

int
SpeciesDataManager::setMotifNetwork(const char* aPtr)
{
	ifstream inFile(aPtr);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1024);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string tfName;
		string tgtName;
		double score=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tfName.append(tok);
			}
			else if(tokCnt==1)
			{
				tgtName.append(tok);
			}	
			else if(tokCnt==2)
			{
				score=atof(tok);				
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		int tfID=varMgr->getVarID(tfName.c_str());
		int tgtID=varMgr->getVarID(tgtName.c_str());
		map<int,double>* tgtSet=NULL;
		if(motifNetwork.find(tfID)==motifNetwork.end())
		{
			tgtSet=new map<int,double>;
			motifNetwork[tfID]=tgtSet;
		}
		else
		{
			tgtSet=motifNetwork[tfID];
		}
		(*tgtSet)[tgtID]=score;
	}
	inFile.close();
	return 0;
}


int 
SpeciesDataManager::readRegulators(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string regulator(buffer);
		regulatorSet[regulator]=0;
	}
	inFile.close();
	return 0;
}

int 
SpeciesDataManager::readTargets(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string target(buffer);
		targetSet[target]=0;
	}
	inFile.close();
	return 0;
}
	
VariableManager*
SpeciesDataManager::getVariableManager()
{
	return varMgr;
}

FactorGraph*
SpeciesDataManager::getFactorGraph()
{
	return fgraph;
}

EvidenceManager*
SpeciesDataManager::getEvidenceManager()
{
	return evMgr;
}

PotentialManager*
SpeciesDataManager::getPotentialManager()
{
	return potMgr;
}

const char*
SpeciesDataManager::getOutputLoc()
{
	return outputLoc;
}

map<string,int>& 
SpeciesDataManager::getRegulators()
{
	return regulatorSet;
}

map<string,int>&
SpeciesDataManager::getTargets()
{
	return targetSet;
}

map<int,map<int,double>*>&
SpeciesDataManager::getMotifNetwork()
{
	return motifNetwork;
}
