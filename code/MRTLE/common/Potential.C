#include <iostream>
#include <math.h>
#include "Evidence.H"

#include "Variable.H"
#include "Potential.H"

#include "gsl/gsl_randist.h"

Potential::Potential()
{
	mean=NULL;
	covariance=NULL;
	inverse=NULL;
	matrixInd=0;
	mbcondVar=0;
	mbcondMean_Part=0;
}

Potential::~Potential()
{
	varSet.clear();
	factorVariables.clear();
	markovBlnktVariables.clear();
	if(covariance!=NULL)
	{
		delete covariance;
	}
	if(inverse!=NULL)
	{
		delete inverse;
	}
	if(mean!=NULL)
	{
		delete mean;
	}
	meanWrap.clear();
	vIDMatIndMap.clear();
	matIndvIDMap.clear();
	mbcondMean_Vect.clear();
}

int
Potential::resetVarSet()
{
	matrixInd=0;
	varSet.clear();
	factorVariables.clear();
	markovBlnktVariables.clear();
	meanWrap.clear();
	vIDMatIndMap.clear();
	matIndvIDMap.clear();
	return 0;
}


int
Potential::initMatrix(int k)
{
	covariance=new Matrix(k,k);
	mean=new Matrix(k,1);
	return 0;
}

int 
Potential::setAssocVariable(Variable* var,Potential::VariableRole vRole)
{
	varSet[var->getID()]=var;
	vIDMatIndMap[var->getID()]=matrixInd;
	matIndvIDMap[matrixInd]=var->getID();
	matrixInd++;
	switch(vRole)
	{
		case Potential::FACTOR:
		{
			factorVariables[var->getID()]=0;
			break;
		}
		case Potential::MARKOV_BNKT:
		{
			markovBlnktVariables[var->getID()]=0;
			break;
		}
	}
	return 0;
}

VSET& 
Potential::getAssocVariables()
{
	return varSet;
}

//The goal of this function is to simply change the variable on which we must condition
//on when computing the conditional potential
int
Potential::updateFactorVariable(int newfVar)
{
	if(factorVariables.size()!=1)
	{
		cout << "Too many factor variables " << endl;
		return -1;
	}
	int currfVar=factorVariables.begin()->first;
	if(currfVar==newfVar)
	{
		return 0;
	}
	Variable* fVar=varSet[currfVar];
	INTINTMAP_ITER vIter=markovBlnktVariables.find(newfVar);
	if(vIter==markovBlnktVariables.end())
	{
		cout <<"Did not find variable " << newfVar <<" in the mb set " << endl;
		return -1;
	}
	markovBlnktVariables.erase(vIter);
	markovBlnktVariables[currfVar]=0;
	factorVariables.clear();
	factorVariables[newfVar]=0;
	return 1;
}

int
Potential::potZeroInit()
{
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		meanWrap[vIter->first]=0;
	}
	int row=varSet.size();
	covariance=new Matrix(row,row);
	covariance->setAllValues(0);
	mean=new Matrix(row,1);
	mean->setAllValues(0);
	return 0;
}


int
Potential::potZeroInit_MeanOnly()
{
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		meanWrap[vIter->first]=0;
	}
	return 0;
}

//Initialize the mean and variance using user specified values
int
Potential::potCustomInit(double mVal, double varVal)
{
	covariance=new Matrix(varSet.size(),varSet.size());
	mean=new Matrix(varSet.size(),1);
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		meanWrap[vIter->first]=mVal;
		int ind=vIDMatIndMap[vIter->first];
		mean->setValue(mVal,ind,0);
	}
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		int i=vIDMatIndMap[vIter->first];
		for(VSET_ITER uIter=varSet.begin();uIter!=varSet.end();uIter++)
		{
			int j=vIDMatIndMap[uIter->first];
			covariance->setValue(varVal,i,j);
		}
	}
	return 0;
}

//Get the joint prob value for a particular configuration
double 
Potential::getJointPotValueFor(INTDBLMAP& varConf)
{
	string aKey;
	double pVal=0;
	Matrix* valMat=new Matrix(varSet.size(),1);
	for(INTDBLMAP_ITER idIter=varConf.begin();idIter!=varConf.end();idIter++)
	{
		int i=vIDMatIndMap[idIter->first];
		valMat->setValue(idIter->second,i,0);
	}
	Matrix* meanDiff=valMat->subtractMatrix(mean);
	Matrix* diffT=meanDiff->transMatrix();
	Matrix* p1=diffT->multiplyMatrix(inverse);
	Matrix* p2=p1->multiplyMatrix(meanDiff);
	double prod=p2->getValue(0,0);
	pVal=exp(-0.5*prod);
	pVal=pVal/normFactor;
	delete meanDiff;
	delete diffT;
	delete p1;
	delete p2;
	return pVal;
}


double 
Potential::getJointPotValueFor(STRINTMAP& varConf)
{
	cout <<"Not implemented " << endl;
	return 0;
}


double 
Potential::getJointPotValueForConf(string& varConf)
{
	cout <<"Not implemented " << endl;
	return 0;
}


int 
Potential::updateMean(int vID,double mVal)
{
	int mID=vIDMatIndMap[vID];
	meanWrap[vID]=mVal;
	mean->setValue(mVal,mID,0);
	return 0;
}

int 
Potential::updateCovariance(int vID,int uID,double sVal)
{
	int mID=vIDMatIndMap[vID];
	int nID=vIDMatIndMap[uID];
	covariance->setValue(sVal,mID,nID);
	return 0;
}


int
Potential::dumpPotential(ostream& oFile)
{
	oFile <<"Mean"<< endl;
	if(meanWrap.size()>0)
	{
		for(INTDBLMAP_ITER sIter=meanWrap.begin();sIter!=meanWrap.end();sIter++)
		{
			oFile <<varSet[sIter->first]->getName() << "\t" << sIter->second << endl;
		}
	}
	return 0;
}

int
Potential::showVariables(ostream& oFile)
{
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		oFile <<" " << vIter->second->getName();
	}
	oFile << endl;
	return 0;
}

//The mean vector is ok. Now we need to get the invariance of the covariance and determinant.
//Also compute the normalization factor
int
Potential::makeValidJPD()
{
	//covariance->showMatrix(cout);
	inverse=covariance->invMatrix();
	//inverse->showMatrix(cout);
	determinant=covariance->detMatrix();
	if(determinant <0)
	{
		cout <<"Negative Determinant " << determinant << endl;
	}
	double n=((double)varSet.size())/2.0;
	normFactor=pow(2*PI,n)*determinant;
	normFactor=sqrt(normFactor);
	return 0;
}


int
Potential::makeValidJPD(gsl_matrix* ludecomp, gsl_permutation* p)
{
	if(inverse!=NULL)
	{
		delete inverse;
	}
	inverse=covariance->invMatrix(ludecomp,p);
	determinant=covariance->detMatrix(ludecomp,p);
	if(determinant <0)
	{
		cout <<"Negative Determinant " << determinant << endl;
		covariance->showMatrix();
		for(INTINTMAP_ITER vIter=matIndvIDMap.begin();vIter!=matIndvIDMap.end();vIter++)
		{
			int vId=vIter->second;
			cout << vIter->first << " " << vId << " "<< varSet[vId]->getName() << endl;
		}
	}
	double n=((double)varSet.size())/2.0;
	normFactor=pow(2*PI,n)*determinant;
	normFactor=sqrt(normFactor);
	return 0;
}

//This is the conditional entropy and is calculated as H(FVar|MarkovBlnkVar)
int
Potential::calculateEntropy()
{
	double entropy=0;
	double jtEntropy=0;
	//Need to create a smaller version of the covariance matrix discarding the row and columns for the
	//factorVariable
	Matrix* neighbourCov=new Matrix(markovBlnktVariables.size(),markovBlnktVariables.size());
	int fmind=vIDMatIndMap[factorVariables.begin()->first];
	for(INTINTMAP_ITER vIter=markovBlnktVariables.begin();vIter!=markovBlnktVariables.end();vIter++)
	{
		int i=vIDMatIndMap[vIter->first];
		for(INTINTMAP_ITER uIter=markovBlnktVariables.begin();uIter!=markovBlnktVariables.end();uIter++)
		{
			int j=vIDMatIndMap[uIter->first];
			double cv=covariance->getValue(i,j);
			if((i==fmind) || (j==fmind))
			{
				cout <<"Factor and Markov Blanket variables are the same!!" << endl;
				return -1;
			}
			if(i>fmind)
			{
				i--;
			}
			if(j>fmind)
			{
				j--;
			}
			neighbourCov->setValue(cv,i,j);
		}
	}
	double nDet=neighbourCov->detMatrix();
	double commFact=1+log(2*PI);
	double p=((double)markovBlnktVariables.size());
	double mbEntropy=0.5*((p*commFact) + log(nDet));
	double n=((double)varSet.size());
	jointEntropy=0.5*((n*commFact) + log(determinant));
	potEntropy=jointEntropy-mbEntropy;
	return 0;
}


int
Potential::calculateJointEntropy()
{
	double commFact=1+log(2*PI);
	double n=((double)varSet.size());
	jointEntropy=0.5*((n*commFact) + log(determinant));
	return 0;
}


double
Potential::getEntropy()
{
	return potEntropy;
}

double
Potential::getJointEntropy()
{
	return jointEntropy;
}


double
Potential::generateSample(INTDBLMAP& jointConf, int vId,gsl_rng* r)
{
	if(jointConf.find(factorVariables.begin()->first)==jointConf.end())
	{
		cout <<"Fatal error! No variable assignment for " << factorVariables.begin()->first << endl;
		exit(0);
	}
	double newmean=0;
	for(INTDBLMAP_ITER aIter=mbcondMean_Vect.begin();aIter!=mbcondMean_Vect.end();aIter++)
	{
		if(jointConf.find(aIter->first)==jointConf.end())
		{
			cout <<"Fatal error! No variable assignment for " << aIter->first << endl;
			exit(0);
		}
		double aval=jointConf[aIter->first];
		newmean=newmean+(aval*aIter->second);
	}
	newmean=newmean+mbcondMean_Part;
	double x=gsl_ran_gaussian(r,sqrt(mbcondVar));
	x=x+newmean;
	return x;
}


int 
Potential::copyMe(Potential** apot)
{
	*apot=new Potential;
	for(INTINTMAP_ITER vIter=matIndvIDMap.begin();vIter!=matIndvIDMap.end();vIter++)
	{
		Variable* v=varSet[vIter->second];
		if(factorVariables.find(vIter->second)!=factorVariables.end())
		{
			(*apot)->setAssocVariable(v,Potential::FACTOR);
		}
		else
		{
			(*apot)->setAssocVariable(v,Potential::MARKOV_BNKT);
		}
	}
	
	for(INTDBLMAP_ITER sdIter=meanWrap.begin();sdIter!=meanWrap.end();sdIter++)
	{
		(*apot)->updateMean(sdIter->first,sdIter->second);
	}
	for(INTDBLMAP_ITER uIter=meanWrap.begin();uIter!=meanWrap.end();uIter++)
	{
		int i=vIDMatIndMap[uIter->first];
		for(INTDBLMAP_ITER vIter=meanWrap.begin();vIter!=meanWrap.end();vIter++)
		{
			int j=vIDMatIndMap[vIter->first];
			double cv=covariance->getValue(i,j);
			(*apot)->updateCovariance(i,j,cv);
		}
	}
	return 0;
}


int 
Potential::initMBCovMean()
{
	int vId=factorVariables.begin()->first;
	int vIdmId=vIDMatIndMap[vId];
	mbcondVar=covariance->getValue(vIdmId,vIdmId);
	mbcondMean_Part=meanWrap[vId];
	if(markovBlnktVariables.size()==0)
	{
		return 0;
	}
	Matrix* mbcov=new Matrix(markovBlnktVariables.size(),markovBlnktVariables.size());
	Matrix* mbmargvar=new Matrix(1,markovBlnktVariables.size());
	INTINTMAP localMatIDMap;
	for(INTDBLMAP_ITER uIter=meanWrap.begin();uIter!=meanWrap.end();uIter++)
	{
		int i=vIDMatIndMap[uIter->first];
		int inew=i;
		if(i>vIdmId)
		{
			inew--;
		}
		for(INTDBLMAP_ITER vIter=meanWrap.begin();vIter!=meanWrap.end();vIter++)
		{
			if(vIter->first==vId)
			{
				continue;
			}
			int j=vIDMatIndMap[vIter->first];
			double cv=covariance->getValue(i,j);
			if(j>vIdmId)
			{
				j--;
			}
			if(uIter->first==vId)
			{
				mbmargvar->setValue(cv,0,j);
			}
			else
			{
				mbcov->setValue(cv,inew,j);
			}
		}
		if(uIter->first!=vId)
		{
			localMatIDMap[inew]=uIter->first;
		}
	}
	Matrix* covInv=mbcov->invMatrix();
	Matrix* prod1=mbmargvar->multiplyMatrix(covInv);
	for(INTINTMAP_ITER aIter=localMatIDMap.begin();aIter!=localMatIDMap.end();aIter++)
	{
		double aVal=prod1->getValue(0,aIter->first);
		double bVal=mbmargvar->getValue(0,aIter->first);
		mbcondVar=mbcondVar-(aVal*bVal);
		mbcondMean_Vect[aIter->second]=aVal;
		double cVal=meanWrap[aIter->second];
		mbcondMean_Part=mbcondMean_Part-(cVal*aVal);
	}
	localMatIDMap.clear();
	if(mbcondVar<1e-10)
	{
		//mbcondVar=1e-10;
	}
	delete mbcov;
	delete mbmargvar;
	delete covInv;
	delete prod1;
	return 0;
}

double 
Potential::getCondPotValueFor(INTDBLMAP& assignment)
{
	if(assignment.find(factorVariables.begin()->first)==assignment.end())
	{
		cout <<"Fatal error! No variable assignment for " << factorVariables.begin()->first << endl;
		exit(0);
	}
	double newmean=0;
	for(INTDBLMAP_ITER aIter=mbcondMean_Vect.begin();aIter!=mbcondMean_Vect.end();aIter++)
	{
		if(assignment.find(aIter->first)==assignment.end())
		{
			cout <<"Fatal error! No variable assignment for " << aIter->first << endl;
			exit(0);
		}
		double aval=assignment[aIter->first];
		newmean=newmean+(aval*aIter->second);
	}
	newmean=newmean+mbcondMean_Part;
	double normsq=2*PI*mbcondVar;
	double norm=sqrt(2*PI*mbcondVar);
	norm=sqrt(normsq);
	double x=assignment[factorVariables.begin()->first];
	double dev=(x-newmean)*(x-newmean);
	dev=dev/(2*mbcondVar);
	double eval=exp(-1.0*dev);
	double pval=eval/norm;
	return pval;
}



double 
Potential::getCondPotValueFor(map<int,Evidence*>* evidMap)
{
	if(evidMap->find(factorVariables.begin()->first)==evidMap->end())
	{
		cout <<"Fatal error! No variable assignment for " << factorVariables.begin()->first << endl;
		exit(0);
	}
	double newmean=0;
	for(INTDBLMAP_ITER aIter=mbcondMean_Vect.begin();aIter!=mbcondMean_Vect.end();aIter++)
	{
		if(evidMap->find(aIter->first)==evidMap->end())
		{
			cout <<"Fatal error! No variable assignment for " << aIter->first << endl;
			exit(0);
		}
		Evidence* evid=(*evidMap)[aIter->first];
		double aval=evid->getEvidVal();
		newmean=newmean+(aval*aIter->second);
	}
	newmean=newmean+mbcondMean_Part;
	double normsq=2*PI*mbcondVar;
	double norm=sqrt(2*PI*mbcondVar);
	norm=sqrt(normsq);
	Evidence* fevid=(*evidMap)[factorVariables.begin()->first];
	double x=fevid->getEvidVal();
	double dev=(x-newmean)*(x-newmean);
	dev=dev/(2*mbcondVar);
	double eval=exp(-1.0*dev);
	double pval=eval/norm;
	return pval;
}

double
Potential::getCondVariance()
{
	return mbcondVar;
}

double
Potential::getCondBias()
{
	return mbcondMean_Part;
}

INTDBLMAP&
Potential::getCondWeight()
{
	return mbcondMean_Vect;
}
