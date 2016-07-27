#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "SpeciesDistance.H"

SpeciesDistance::SpeciesDistance()
{
	root=NULL;
}

SpeciesDistance::~SpeciesDistance()
{
}

int 
SpeciesDistance::readSpeciesTree(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	//double p_gain = pgain;
	//double p_maintain_edge = pmaintain;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string childSpeciesName;
		string parentSpeciesName;
		string childtype;
		double p_gain=0;
		double p_maintain_edge=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				childSpeciesName.append(tok);	
			}
			else if(tokCnt==1)
			{
				childtype.append(tok);
			}
			else if(tokCnt==2)
			{
				parentSpeciesName.append(tok);
			}
			else if(tokCnt==3)
			{
				p_gain=atof(tok);
			}
			else if(tokCnt==4)
			{
				p_maintain_edge=1-atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		SpeciesDistance::Species* childSpecies=NULL;
		SpeciesDistance::Species* parentSpecies=NULL;
		if(speciesSet.find(childSpeciesName)==speciesSet.end())
		{
			childSpecies=new SpeciesDistance::Species;
			childSpecies->name.append(childSpeciesName.c_str());
			childSpecies->parent=NULL;
			childSpecies->leftchild=NULL;
			childSpecies->rightchild=NULL;
			speciesSet[childSpeciesName]=childSpecies;
		}
		else
		{
			childSpecies=speciesSet[childSpeciesName];
		}
		if(speciesSet.find(parentSpeciesName)==speciesSet.end())
		{
			parentSpecies=new SpeciesDistance::Species;
			parentSpecies->name.append(parentSpeciesName.c_str());
			parentSpecies->parent=NULL;
			parentSpecies->leftchild=NULL;
			parentSpecies->rightchild=NULL;
			speciesSet[parentSpeciesName]=parentSpecies;
		}
		else
		{
			parentSpecies=speciesSet[parentSpeciesName];
		}
		childSpecies->parent=parentSpecies;
		if(strcmp(childtype.c_str(),"left")==0)
		{
			parentSpecies->leftchild=childSpecies;
		}	
		else if(strcmp(childtype.c_str(),"right")==0)
		{
			parentSpecies->rightchild=childSpecies;
		}
		else
		{
			cout <<"Unknown edge type " << childtype.c_str() << endl;
		}
		if(parentSpecies->parent==NULL)
		{
			root=parentSpecies;
		}
		//Probabilities are for the child node
		childSpecies->p_maintain_edge=p_maintain_edge;
		childSpecies->p_gain=p_gain;
	}

	inFile.close();
	return 0;
}


//The probability that an edge is maintained in a child given that the edge is present in the ancestor
double 
SpeciesDistance::getProbGainGain(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	double gain_gain=species->p_maintain_edge;
	return gain_gain;
}

//The probability that an edge is gained in a child given that the edge is absent in the ancestor
double 
SpeciesDistance::getProbGainLoss(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	double gain_loss=species->p_gain;
	return gain_loss;
}

//The probability than an edge is lost in a child given that the edge was present in the ancestor
double 
SpeciesDistance::getProbLossGain(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	double loss_gain=1-species->p_maintain_edge;
	return loss_gain;
}

//The  probability that an edge is lost in a child given the edge was absent in the ancestor
double 
SpeciesDistance::getProbLossLoss(string& spName)
{
	if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}
	Species* species=speciesSet[spName];
	double loss_loss=1-species->p_gain;
	return loss_loss;
}

double
SpeciesDistance::getEdgeStatusProb(map<string,int>& edgeStatus)
{
	/*for(map<string,int>::iterator eIter=edgeStatus.begin();eIter!=edgeStatus.end();eIter++)
	{
		cout <<eIter->first <<"\t" << eIter->second << endl;
	}*/
	double edgeprior=0;
	int binkey=0;
	int binpos=0;
	for(map<string,int>::iterator eIter=edgeStatus.begin();eIter!=edgeStatus.end();eIter++)
	{
		if(eIter->second!=0)
		{	
			int val=(int)pow(2,binpos);
			binkey=binkey+val;
		}
		binpos=binpos+1;
	}
	if(binEdgeKeyProbMap.find(binkey)!=binEdgeKeyProbMap.end())
	{
		return binEdgeKeyProbMap[binkey];
	}
	//Start with the root
	bool edgeInParent[]={false,true};
	
	for(int i=0;i<2;i++)
	{
		double leftScore=getSubTree(edgeInParent[i],root->leftchild,edgeStatus);
		double rightScore=getSubTree(edgeInParent[i],root->rightchild,edgeStatus);
		//changed 6/4 to test setting the prior probability of an edge in the root node = .2
		double score;
		if(i==0)
		{
			score = 0.8*leftScore*rightScore;
		}
		else
		{
			score = 0.2*leftScore*rightScore;
		}		
		edgeprior=edgeprior+score;
	}
	binEdgeKeyProbMap[binkey]=edgeprior;
	return edgeprior;
}

double
SpeciesDistance::getSubTree(bool parentEdge, Species* child, map<string,int>& edgeStatus)
{
	double score=0;
	if(child->leftchild==NULL && child->rightchild==NULL)
	{
		//This is a leaf node
		int status=edgeStatus[child->name];
		if(parentEdge)
		{
			if(status==0)
			{
				score=getProbLossGain(child->name);
			}
			else
			{
				score=getProbGainGain(child->name);
			}
		}
		else
		{	
			if(status==0)
			{
				score=getProbLossLoss(child->name);
			}
			else
			{
				score=getProbGainLoss(child->name);
			}
		}
	}
	else
	{
		//This is an interior node
		bool edgeInParent[]={false,true};
		for(int i=0;i<2;i++)
		{
			double leftScore=getSubTree(edgeInParent[i],child->leftchild,edgeStatus);
			double rightScore=getSubTree(edgeInParent[i],child->rightchild,edgeStatus);
			double prob=0;
			if(parentEdge)
			{
				if(edgeInParent[i])
				{
					prob=getProbGainGain(child->name);
				}
				else
				{
					prob=getProbLossGain(child->name);
				}	
			}
			else
			{
				if(edgeInParent[i])
				{
					prob=getProbGainLoss(child->name);
				}
				else
				{
					prob=getProbLossLoss(child->name);
				}	
				
			}
			score=score+(prob*leftScore*rightScore);
		}
	}
	return score;
}

SpeciesDistance::Species*
SpeciesDistance::getRoot()
{
	return root;
}
