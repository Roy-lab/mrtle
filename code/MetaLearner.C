#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

#include "FactorGraph.H"
#include "MetaMove.H"

#include "Utils.H"

#include "SpeciesDistance.H"
#include "SpeciesDataManager.H"
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "MetaLearner.H"
#include "Matrix.H"
#include <chrono>
#include <unistd.h>
using namespace std::chrono;
using namespace std;

MetaLearner::MetaLearner()
{
    restrictedFName[0]='\0';
    lambdaLearnSize=0;
    learnLambda=false;
    dataPool=false;
    usePLL=true;
    useEntropy=true;
    trueGraphFName[0]='\0';
    preRandomizeSplit=false;
    preRandSeed=-1;
    beta1=-0.9;
    beta2=4.0;
    convThreshold=1e-4;
}

MetaLearner::~MetaLearner()
{
}

int
MetaLearner::setInputFName(const char* aFName)
{
    //strcpy(inputFName,aFName);
    inputFName = aFName;
    return 0;
}

int
MetaLearner::setMaxFactorSize(int aVal)
{
    maxFactorSize=aVal;
    return 0;
}

int
MetaLearner::setMaxFactorSize_Approx(int aVal)
{
    maxFactorSizeApprox=aVal;
    return 0;
}

int 
MetaLearner::setMBCntForApproxGreedy(int aCnt)
{
    maxMBCntGreedy=aCnt;
    return 0;
}

int 
MetaLearner::setPenalty(double aVal)
{
    penalty=aVal;
    return 0;
}

int 
MetaLearner::setConvergenceThreshold(double aVal)
{
    convThreshold=aVal;
    return 0;
}

int
MetaLearner::setRestrictedList(const char* aFName)
{
    strcpy(restrictedFName,aFName);
    ifstream inFile(restrictedFName);
    string buffer;
    while(inFile.good())
    {
        getline(inFile,buffer);
        if(buffer.length()<=0)
        {
            continue;
        }
        int ogid=atoi(buffer.c_str());
        inputRegulatorOGs[ogid]=0;
    }
    inFile.close();
    return 0;
}

int
MetaLearner::setLambdaLearnSize(int vSize)
{
    lambdaLearnSize=vSize;
    return 0;
}

int
MetaLearner::setLambdaLearn()
{
    learnLambda=true;
    return 0;
}

int
MetaLearner::setTrueGraph(const char* aGraphName)
{
    strcpy(trueGraphFName,aGraphName);
    return 0;
}


int 
MetaLearner::setPreRandomizeSplit()
{
    preRandomizeSplit=true;
    return 0;
}

int
MetaLearner::setPreRandomizeSplitSeed(int seed)
{
    preRandSeed=seed;
    return 0;
}


int
MetaLearner::setScoreOpt_Entropy()
{
    usePLL=false;
    useEntropy=true;
    return 0;
}

int
MetaLearner::setScoreOpt_PLL()
{
    usePLL=true;
    useEntropy=false;
    return 0;
}


int 
MetaLearner::setInputOGList(const char* aFName)
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
        int og=atoi(buffer);
        inputOGList[og]=0;
    }
    inFile.close();
    return 0;
}

int 
MetaLearner::setOrthogroupReader(MappedOrthogroupReader* aPtr)
{
    ogr=aPtr;
    return 0;
}

int 
MetaLearner::setSpeciesDistances(SpeciesDistance* aPtr)
{
    speciesData=aPtr;
    return 0;
}

int
MetaLearner::setBeta1(double b1)
{
    beta1 = b1;
    return 0;
}

int
MetaLearner::setBeta2(double b2)
{
    beta2 = b2;
    return 0;
}

int
MetaLearner::initOld()
{
    /*
    ifstream inFile(inputFName);
    char buffer[1024];
    int datasetId=1;
    while(inFile.good())
    {
        inFile.getline(buffer,1023);
        
        if(strlen(buffer)<=0)
        {
            continue;
        }
        if(strstr(buffer,"#")!=NULL)
        {
            continue;
        }
        int tokCnt=0;
        string specName;
        char datasetSuff[1024];
        char outputLoc[1024];
        char regulatorFName[1024];
        char targetFName[1024];
        char motifNetwork[1024];
        char* tok=strtok(buffer,"\t");
        while(tok!=NULL)
        {
            if(tokCnt==0)
            {
                specName.append(tok);
            }
            else if(tokCnt==1)
            {
                strcpy(datasetSuff,tok);
            }
            else if(tokCnt==2)
            {
                strcpy(outputLoc,tok);
            }
            else if(tokCnt==3)
            {
                strcpy(regulatorFName,tok);
            }
            else if(tokCnt==4)
            {
                strcpy(targetFName,tok);
            }
            else if(tokCnt==5)
            {
                strcpy(motifNetwork,tok);   
            }
            tok=strtok(NULL,"\t");
            tokCnt++;
        }

        PotentialManager* potMgr=new PotentialManager;
        EvidenceManager* evMgr=new EvidenceManager;
        if(preRandomizeSplit)
        {
            evMgr->setPreRandomizeSplit();
            evMgr->setPreRandomizeSplitSeed(preRandSeed);
        }

        potMgr->setEvidenceManager(evMgr);
        
        char fName[256];
        sprintf(fName,"%s.model",datasetSuff);
        VariableManager* varMgr=new VariableManager;
        Error::ErrorCode eCode=varMgr->readVariables(fName);
        if(eCode!=Error::SUCCESS)
        {
            cout << Error::getErrorString(eCode) << endl;
            return -1;
        }
        evMgr->setVariableManager(varMgr);
        sprintf(fName,"%s.data",datasetSuff);
        eCode=evMgr->loadEvidenceFromFile_Continuous(fName);
        if(eCode!=Error::SUCCESS)
        {
            cout <<Error::getErrorString(eCode) << endl;
            return -1;
        }
        potMgr->setOutputDir(outputLoc);
        SpeciesDataManager* spMgr=new SpeciesDataManager;
        speciesDataSet[specName]=spMgr;
        spMgr->setVariableManager(varMgr);
        spMgr->setEvidenceManager(evMgr);
        spMgr->setPotentialManager(potMgr);
        spMgr->createFactorGraph();
        spMgr->setOutputLoc(outputLoc);
        spMgr->readRegulators(regulatorFName);
        spMgr->readTargets(targetFName);
        spMgr->setMotifNetwork(motifNetwork);
        speciesIDNameMap[datasetId]=specName;
        speciesNameIDMap[specName]=datasetId;
        speciesDataSet[specName]=spMgr;
        datasetId=datasetId*2;
    }
    inFile.close();
    cout <<"Read data for " << speciesDataSet.size() << " datasets" << endl;
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        cout << sIter->first << endl;
    }
    */
    return 0;
}

int
MetaLearner::init()
{
    ifstream inFile(inputFName);
    string buffer;
    //char buffer[1024];
    int datasetId = 1;  
    while (inFile.good())
    {
        getline(inFile, buffer);
        if (buffer.empty() || buffer.find("#") == 0)
        {
            continue;
        }
        /*
        inFile.getline(buffer,1023);        
        if(strlen(buffer)<=0)
        {
            continue;
        }
        if(strstr(buffer,"#")!=NULL)
        {
            continue;
        }
        */
        vector<string> strs;

        // strip trailing newlines
        buffer.erase(buffer.find_last_not_of(" \n\r\t") + 1);

        strs = Utils::split(buffer, '\t');
        assert(strs.size() == 6);
        string specName = strs[0];
        string datasetSuff = strs[1];
        string outputLoc = strs[2];
        string regulatorFName = strs[3];
        string targetFName = strs[4];
        string motifNetwork = strs[5];

        PotentialManager* potMgr = new PotentialManager;
        EvidenceManager* evMgr = new EvidenceManager;
        if (preRandomizeSplit)
        {
            evMgr->setPreRandomizeSplit();
            evMgr->setPreRandomizeSplitSeed(preRandSeed);
        }

        potMgr->setEvidenceManager(evMgr);
        
        string tableFileName = datasetSuff ;//+ ".table";
        readEvidenceTable(tableFileName);

        /*VariableManager* varMgr = new VariableManager;
        Error::ErrorCode eCode = varMgr->readVariablesFromTable(inputTable);
        if(eCode != Error::SUCCESS)
        {
            cout << Error::getErrorString(eCode) << endl;
            return -1;
        }
        evMgr->setVariableManager(varMgr);*/
        Error::ErrorCode eCode=evMgr->loadEvidenceFromTable(inputTable);
        if(eCode!=Error::SUCCESS)
        {
            cout <<Error::getErrorString(eCode) << endl;
            return -1;
        }
        VariableManager* varMgr=evMgr->getVariableManager();
        potMgr->setOutputDir(outputLoc.c_str());
        SpeciesDataManager* spMgr=new SpeciesDataManager;
        speciesDataSet[specName]=spMgr;
        spMgr->setVariableManager(varMgr);
        spMgr->setEvidenceManager(evMgr);
        spMgr->setPotentialManager(potMgr);
        spMgr->createFactorGraph();
        spMgr->setOutputLoc(outputLoc.c_str());
        spMgr->readRegulators(regulatorFName.c_str());
        spMgr->readTargets(targetFName.c_str());
        spMgr->setMotifNetwork(motifNetwork.c_str());
        speciesIDNameMap[datasetId]=specName;
        speciesNameIDMap[specName]=datasetId;
        speciesDataSet[specName]=spMgr;
        cout << specName << ": " << datasetId << endl;
        datasetId=datasetId*2;
    }
    inFile.close();
    cout <<"Read data for " << speciesDataSet.size() << " datasets" << endl;
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        cout << sIter->first << endl;
    }
    return 0;
}

void
MetaLearner::readEvidenceTable(string fileName)
{
    /* Reads gene names and expression levels from a tab-separated file, where the first row is assumed (for now)
       to be headers. In subsequent rows, the first column is a gene name and remaining columns are expression levels.
       This method just loads the meaningful lines into a vector that is then passed to VariableManager::readVariablesFromTable()
       and EvidenceManager::loadEvidenceFromTable(). It's done this way so as not to break encapsulation of private members
       of VariableManager and EvidenceManager.
    */
    cout << "readEvidenceTable:" << fileName << endl;
    ifstream inFile(fileName);
    assert(inFile);

    bool headerLine = true;
    string inputLine;
    inputTable.clear();

    // Read all the lines with meaningful data
    while (getline(inFile, inputLine))
    {
        if (inputLine.empty() || inputLine.find('#') == 0)
        {
            continue;
        }
        if (headerLine)
        {
            headerLine = false;
            continue;
        }
        inputTable.push_back(inputLine);
    }
}
//shilu: train on the whole dataset
int
MetaLearner::doCrossValidation()
{
    totalFoldCnt=1;
    /*gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
    for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    {
        EvidenceManager* evMgr=eIter->second->getEvidenceManager();
        evMgr->setFoldCnt(foldCnt);
        //evMgr->splitData(0);
    }
    gsl_rng_free(r);*/
    //The first key is for the fold number
    //For each fold we have a trained model. For each trained model we have the likelihood on 
    //all the test sets, including the self test.
    int f=0;
    for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    {
        EvidenceManager* evMgr=eIter->second->getEvidenceManager();
        evMgr->setFoldCnt(totalFoldCnt);
        evMgr->splitData(f);
        PotentialManager* potMgr=eIter->second->getPotentialManager();
        potMgr->reset();
        potMgr->init(); //potMgr->init(f);
        char outputDir[1024];
        sprintf(outputDir,"%s/fold%d",eIter->second->getOutputLoc(),f);
        char foldOutputDirCmd[1024];
        sprintf(foldOutputDirCmd,"mkdir -p %s",outputDir);
        system(foldOutputDirCmd);
    }
    //clearFoldSpecData();
    start(f);
    char scoreFName[1024];
    sprintf(scoreFName,"%s/scoreFile.txt",speciesDataSet.begin()->second->getOutputLoc());
    ofstream sFile(scoreFName);
    for(map<int,double>::iterator dIter=finalScores.begin();dIter!=finalScores.end();dIter++)
    {
        sFile<< dIter->second << endl;
    }
    sFile.close();
    return 0;
}

int
MetaLearner::doCrossValidation(int foldCnt)
{
    totalFoldCnt=foldCnt;
    gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
    for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    {
        EvidenceManager* evMgr=eIter->second->getEvidenceManager();
        evMgr->setFoldCnt(foldCnt);
        //evMgr->splitData(0);
    }
    gsl_rng_free(r);
    //The first key is for the fold number
    //For each fold we have a trained model. For each trained model we have the likelihood on 
    //all the test sets, including the self test.
    for(int f=0;f<foldCnt;f++)
    {   
        for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
        {
            EvidenceManager* evMgr=eIter->second->getEvidenceManager();
            evMgr->splitData(f);
            PotentialManager* potMgr=eIter->second->getPotentialManager();
            potMgr->reset();
            potMgr->init(f);
            char outputDir[1024];
            sprintf(outputDir,"%s/fold%d",eIter->second->getOutputLoc(),f);
            char foldOutputDirCmd[1024];
            sprintf(foldOutputDirCmd,"mkdir -p %s",outputDir);
            system(foldOutputDirCmd);
        }
        //clearFoldSpecData();
        start(f);
        //Now get all the PLL scores on the validation sets
        for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
        {   
            double pll=0;
            VariableManager* varMgr=eIter->second->getVariableManager();
            VSET& varSet=varMgr->getVariableSet();
            int specID=speciesNameIDMap[eIter->first];
            for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
            {
                pll=pll+getValidationPLLScore_Condition(specID,vIter->first);
            }
            INTDBLMAP* pllSet=NULL;
            if(validationPLLs.find(specID)==validationPLLs.end())
            {
                pllSet=new INTDBLMAP;
                validationPLLs[specID]=pllSet;
            }
            else
            {
                pllSet=validationPLLs[specID];
            }
            (*pllSet)[f]=pll;
        }
    }
    cout <<"Validation PLLs" << endl;
    for(map<int,INTDBLMAP*>::iterator eIter=validationPLLs.begin();eIter!=validationPLLs.end();eIter++)
    {
        INTDBLMAP* pllSet=eIter->second;
        cout <<"Vset"<<eIter->first;
        for(INTDBLMAP_ITER fIter=pllSet->begin();fIter!=pllSet->end();fIter++)
        {
            cout <<"\t" << fIter->second;
        }
        cout << endl;
    }
    char scoreFName[1024];
    sprintf(scoreFName,"%s/scoreFile.txt",speciesDataSet.begin()->second->getOutputLoc());
    ofstream sFile(scoreFName);
    for(map<int,double>::iterator dIter=finalScores.begin();dIter!=finalScores.end();dIter++)
    {
        sFile<< dIter->second << endl;
    }
    sFile.close();
    return 0;
}

int
MetaLearner::start(int f)
{
    cout << "MetaLearner::start beta1=" << beta1 << " beta2=" << beta2 << endl;
    auto start = high_resolution_clock::now();
    //Repeat until convergence
    double currGlobalScore=-1;
    int maxMBSizeApprox=maxFactorSizeApprox-1;
    //int currK=1;
    int currK=maxMBSizeApprox;
    double initPrior=precomputeEmptyGraphPrior();
    precomputePerSpeciesPrior();
    initEdgeSet();
    if(strlen(trueGraphFName)==0)
    {
        while(currK<=maxMBSizeApprox)
        {
            int iter = 0; // vperiyasamy added
            bool notConverged=true;
            while(notConverged && iter < 50) // hardcode 50 iterations for now
            {
                cout << "ITERATION " << iter << endl;
                //collect the candidate edges
                collectMoves_Orthogroups(currK);
                //sortMoves();
                double diff0=makeMoves();
                if(currGlobalScore==-1)
                {
                    currGlobalScore=getInitScore();
                //  currGlobalScore=currGlobalScore+initPrior;
                }
                double priorChange=getPriorDelta();
                //collectMoves_Deletions();
                //sortMoves();
                //makeDelMoves();
                //priorChange=getPriorDelta();
                double newScore=getScore();
                cout <<"newScore=" << newScore << " diff+currGlobalScore=" << diff0+currGlobalScore << " difference=" << newScore-diff0-currGlobalScore << endl;
                //newScore=newScore+initPrior+priorChange;
                double diff=0;
                diff=newScore-currGlobalScore;
                if(diff<=convThreshold)
                {
                    notConverged=false;
                }
                //dumpAllGraphs(currK,f); //comment out by shilu, no need to output for each iteration
                currGlobalScore=newScore;
                cout << "ITERATION " << iter << " Current "<< currK <<" Score=" << currGlobalScore << " diff=" << diff<< endl;
                iter++;
            }
            dumpAllGraphs(currK,f);
            currK++;
        }
        cout <<"Final Score " << currGlobalScore << endl;
        finalScores[f]=currGlobalScore;
        showModelParameters(f);
    }
    else
    {
        //populateGraphsFromFile();
        double globalScore=getScore();
        generateData(10000,30000,f);
        showModelParameters(f);
        //cout <<"Final Score " << globalScore<< endl;
        cout <<"Final Score\tTOTALFOLD" << totalFoldCnt <<"_fold"<< f<<"\t" << globalScore<< endl;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);    
    cout << "MetaLearner::start runtime: " <<duration.count() << " ms" <<endl;
    return 0;
}


double
MetaLearner::getInitScore()
{
    double initScore=0;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {
        FactorGraph* fg=gIter->second->getFactorGraph();
        map<int,SlimFactor*>& factorSet=fg->getAllFactors();
        for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
        {   
            SlimFactor* sFactor=fIter->second;
            initScore=initScore+sFactor->marginalLL;
        //  cout << gIter->first << "\t" << fIter->first <<"\t" << sFactor->marginalLL<< endl;
        }
    }
    return initScore;
}

double
MetaLearner::getScore()
{
    double gScore=0;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {   
        FactorGraph* fg=gIter->second->getFactorGraph();
        map<int,SlimFactor*>& factorSet=fg->getAllFactors();
        for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
        {
            SlimFactor* sFactor=fIter->second;
            gScore=gScore+sFactor->mbScore;
        }
    }
    return gScore;
}

double
MetaLearner::getPriorDelta()
{
    double oldStructPrior=0;
    double newStructPrior=0;
    //Need to consider the old contribution of the edges, delete that from the overall prior and add the new contribution
    for(map<string,STRINTMAP*>::iterator edgeIter=affectedOGPairs.begin();edgeIter!=affectedOGPairs.end();edgeIter++)
    {
        double aval=speciesData->getEdgeStatusProb(*(edgeIter->second));
        double edgePrior=log(aval);
        double oldEdgePrior=ogpairPrior[edgeIter->first];
        oldStructPrior=oldStructPrior+oldEdgePrior;
        newStructPrior=newStructPrior+edgePrior;
        INTINTMAP* currEdgeStatus=edgeConditionMap[edgeIter->first];
        STRINTMAP* newEdgeStatus=edgeIter->second;
        for(STRINTMAP_ITER sIter=newEdgeStatus->begin();sIter!=newEdgeStatus->end();sIter++)
        {
            int specID=speciesNameIDMap[sIter->first];
            (*currEdgeStatus)[specID]=sIter->second;
        }
        ogpairPrior[edgeIter->first]=edgePrior;
    }
    for(map<string,STRINTMAP*>::iterator edgeIter=affectedOGPairs.begin();edgeIter!=affectedOGPairs.end();edgeIter++)
    {
        edgeIter->second->clear();
        delete edgeIter->second;
    }
    affectedOGPairs.clear();
    double priorDelta=newStructPrior-oldStructPrior;
    return priorDelta;
}

int
MetaLearner::clearFoldSpecData()
{
    //Clear existing graphs
    for(map<string,SpeciesDataManager*>::iterator fIter=speciesDataSet.begin();fIter!=speciesDataSet.end();fIter++)
    {
        delete fIter->second;
    }
    speciesDataSet.clear();
    return 0;
}

int
MetaLearner::initEdgeSet()
{
    if(condsetMap.size()==0)
    {
        if(dataPool)
        {
            initCondsetMap();
        }
        else
        {
            initCondsetMap_Nopool();
        }
    }
    initCondsetMap_Tree(speciesData->getRoot());
    for(map<int,INTINTMAP*>::iterator setIter=condsetMap_Tree.begin();setIter!=condsetMap_Tree.end();setIter++)
    {
        INTINTMAP* cset=setIter->second;
        cout << setIter->first;
        for(INTINTMAP_ITER cIter=cset->begin();cIter!=cset->end();cIter++)
        {
            cout << " "<<cIter->first<<"=" << cIter->second;
        }
        cout << endl;
    }

    for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    {
        FactorGraph* condspecGraph=eIter->second->getFactorGraph();
        int specID=speciesNameIDMap[eIter->first];
        INTINTMAP* cset=getConditionSet(specID);
        string condKey;
        genCondSetKey(*cset,condKey);
        map<int,double>* varNeighborhoodPrior=varNeighborhoodPrior_PerSpecies[eIter->first];
        for(int f=0;f<condspecGraph->getFactorCnt();f++)
        {
            SlimFactor* sFactor=condspecGraph->getFactorAt(f);
            /*if(sFactor->fId==111)
            {
                cout <<"found the target of interest" << endl;
            }*/
            double pll=getPLLScore_Condition((string&)eIter->first,sFactor);
            double priorScore=(*varNeighborhoodPrior)[sFactor->fId];
            sFactor->mbScore=pll+priorScore;
            sFactor->marginalLL=pll+priorScore;
	    //cout << eIter->first  <<" f=" <<f <<" pll=" << pll << " priorScore=" << priorScore << endl;
        }
        PotentialManager* potMgr=eIter->second->getPotentialManager();
        EvidenceManager* evMgr=eIter->second->getEvidenceManager();
        //delete data to save memory:
        if(totalFoldCnt==1){
            cout << "delete evMgr" << endl;
            evMgr->deleteData();
        }
    }
    
    return 0;
}

double
MetaLearner::precomputeEmptyGraphPrior()
{
    map<string,int> ogpairEdgeCntMap;
    map<int,int> ogMaxSpecCnt_TF;
    map<int,int> ogMaxSpecCnt_Tgt;
    STRINTMAP edgeStatus;
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        SpeciesDataManager* sdm=sIter->second;
        VariableManager* vMgr=sdm->getVariableManager();
        VSET& varSet=vMgr->getVariableSet();
        map<string,int>& regulators=sdm->getRegulators();
        map<string,int>& targets=sdm->getTargets();
        for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
        {
            Variable* v=vIter->second;
            if(targets.size()>0 && (targets.find(v->getName())==targets.end()))
            {
                continue;
            }
            int ogno_tgt=ogr->getMappedOrthogroupID(v->getName().c_str(),sIter->first.c_str());
            if(inputOGList.find(ogno_tgt) == inputOGList.end())
            {
                continue;
            }
            if(ogno_tgt==-1)
            {
                cout <<"No Orthogroup for " << v->getName() << " species " << sIter->first << endl;
                continue;
            }
            MappedOrthogroup* ogrp_Tgt=ogr->getMappedOrthogroup(v->getName().c_str(),sIter->first.c_str());
            if(ogrp_Tgt==NULL)
            {
                cout <<"No target by name " << v->getName() << " in " << sIter->first << endl;
            }
            GeneMap* gmap=ogrp_Tgt->getSpeciesHits(sIter->first.c_str());
            map<string,map<string,STRINTMAP*>*>& targetsInSpec=gmap->getGeneSet();
            int targetsWithVals=0;
            for(map<string,map<string,STRINTMAP*>*>::iterator tIter=targetsInSpec.begin();tIter!=targetsInSpec.end();tIter++)
            {
                if(vMgr->getVarID(tIter->first.c_str())==-1)
                {
                    continue;
                }
                targetsWithVals++;
            }
            if(ogMaxSpecCnt_Tgt.find(ogno_tgt)==ogMaxSpecCnt_Tgt.end())
            {
                ogMaxSpecCnt_Tgt[ogno_tgt]=targetsWithVals;
            }
            else
            {
                int currval=ogMaxSpecCnt_Tgt[ogno_tgt];
                if(currval<targetsWithVals)
                {
                    ogMaxSpecCnt_Tgt[ogno_tgt]=targetsWithVals;
                }
            }
            //Now we need to figure out how many edges this ogno_tgt-ogno_tf pair is going to contribute to
            //This is essentially the max number of active targets
            for(map<string,int>::iterator rIter=regulators.begin();rIter!=regulators.end();rIter++)
            {
                int regID=vMgr->getVarID(rIter->first.c_str());
                if(regID==-1)
                {
                    continue;
                }
                int ogno_tf=ogr->getMappedOrthogroupID(rIter->first.c_str(),sIter->first.c_str());
                if(ogno_tf==-1)
                {
                    cout <<"No Orthogroup for " <<rIter->first << " species " << sIter->first << endl;
                    continue;
                }
                MappedOrthogroup* ogrp_TF=ogr->getMappedOrthogroup(rIter->first.c_str(),sIter->first.c_str());
                if(ogrp_TF==NULL)
                {
                    cout <<"No gene " << rIter->first << " in " << sIter->first << endl;
                    continue;
                }
                char key[1024];
                sprintf(key,"%d-%d",ogno_tf,ogno_tgt);
                string keystr(key);
                if(edgeConditionMap.find(keystr)==edgeConditionMap.end())
                {
                    INTINTMAP* edgeSpecAssign=new INTINTMAP;
                    for(map<string,int>::iterator tIter=speciesNameIDMap.begin();tIter!=speciesNameIDMap.end();tIter++)
                    {
                        (*edgeSpecAssign)[tIter->second]=0;
                    }
                    edgeConditionMap[keystr]=edgeSpecAssign;
                }
                GeneMap* tfgenemap=ogrp_TF->getSpeciesHits(sIter->first.c_str());
                map<string,map<string,STRINTMAP*>*>& tfsInSpec=tfgenemap->getGeneSet();
                int tfsWithVals=0;
                for(map<string,map<string,STRINTMAP*>*>::iterator tIter=tfsInSpec.begin();tIter!=tfsInSpec.end();tIter++)
                {
                    if(vMgr->getVarID(tIter->first.c_str())==-1)
                    {
                        continue;
                    }
                    tfsWithVals++;
                }
                if(ogMaxSpecCnt_TF.find(ogno_tf)==ogMaxSpecCnt_TF.end())
                {
                    ogMaxSpecCnt_TF[ogno_tf]=tfsWithVals;
                }
                else
                {
                    int currval=ogMaxSpecCnt_TF[ogno_tf];
                    if(currval<tfsWithVals)
                    {
                        ogMaxSpecCnt_TF[ogno_tf]=tfsWithVals;
                    }
                }
            }
            
        }
        edgeStatus[sIter->first]=0;
    }
    double prior=speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior=log(prior);
    int edgeCnt=0;
    for(map<int,int>::iterator tfIter=ogMaxSpecCnt_TF.begin();tfIter!=ogMaxSpecCnt_TF.end();tfIter++)
    {
        for(map<int,int>::iterator tgtIter=ogMaxSpecCnt_Tgt.begin();tgtIter!=ogMaxSpecCnt_Tgt.end();tgtIter++)
        {
            edgeCnt=edgeCnt+(tfIter->second*tgtIter->second);   
            char key[1024];
            sprintf(key,"%d-%d",tfIter->first,tgtIter->first);
            string keystr(key);
            ogpairPrior[keystr]=logPrior;
        }
    }
    double emptyGraphPrior=edgeCnt*logPrior;
    return emptyGraphPrior;
}

int
MetaLearner::precomputePerSpeciesPrior()
{
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        SpeciesDataManager* sdm=sIter->second;
        VariableManager* vMgr=sdm->getVariableManager();
        VSET& varSet=vMgr->getVariableSet();
        map<string,int>& regulators=sdm->getRegulators();
        map<string,int>& targets=sdm->getTargets();
        map<string,double>* edgePresenceProb=new map<string,double>;
        edgePresenceProb_PerSpecies[sIter->first]=edgePresenceProb;
        map<int,double>* varNeighborhoodPrior=new map<int,double>;
        varNeighborhoodPrior_PerSpecies[sIter->first]=varNeighborhoodPrior;

        for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
        {
            Variable* v=vIter->second;
            if(targets.size()>0 && (targets.find(v->getName())==targets.end()))
            {
                continue;
            }
            int ogno_tgt=ogr->getMappedOrthogroupID(v->getName().c_str(),sIter->first.c_str());
            if(inputOGList.find(ogno_tgt)==inputOGList.end())
            {
                continue;
            }
            //Now we need to figure out how many edges this ogno_tgt-ogno_tf pair is going to contribute to
            //This is essentially the max number of active targets
            for(map<string,int>::iterator rIter=regulators.begin();rIter!=regulators.end();rIter++)
            {
                int regID=vMgr->getVarID(rIter->first.c_str());
                if(regID==-1)
                {
                    continue;
                }
                string edgeKey; 
                edgeKey.append(rIter->first.c_str());
                edgeKey.append("\t");
                edgeKey.append(v->getName().c_str());
                double initPrior=getEdgePrior_PerSpecies(regID,vIter->first,sdm);
                (*edgePresenceProb)[edgeKey]=initPrior;
                if(varNeighborhoodPrior->find(vIter->first)==varNeighborhoodPrior->end())
                {
                    (*varNeighborhoodPrior)[vIter->first]=log(1-initPrior);
                }
                else
                {
                    (*varNeighborhoodPrior)[vIter->first]=(*varNeighborhoodPrior)[vIter->first]+log(1-initPrior);
                }
            }
        }
    }
    
    return 0;
}

int
MetaLearner::initCondsetMap()
{
    int lastcondind=(int)pow(2,speciesDataSet.size());
    for(int i=1;i<lastcondind;i++)
    {
        INTINTMAP* cset=new INTINTMAP;
        int currind=i;
        for(map<int,string>::iterator eIter=speciesIDNameMap.begin();eIter!=speciesIDNameMap.end();eIter++)
        {
            if(currind==0)
            {
                (*cset)[eIter->first]=0;
                continue;
            }
            if((currind%2)==1)
            {
                (*cset)[eIter->first]=1;
            }
            else
            {
                (*cset)[eIter->first]=0;
            }
            currind=currind/2;
        }
        condsetMap[i]=cset;
        string condKey;
        genCondSetKey(*cset,condKey);
        condsetKeyIDMap[condKey]=i;
        condsetIDKeyMap[i]=condKey;
    }
    return 0;
}

int
MetaLearner::initCondsetMap_Tree(SpeciesDistance::Species* node)
{
    if(node->leftchild==NULL && node->rightchild==NULL)
    {
        map<int,int>* cset=new map<int,int>;
        int specID=speciesNameIDMap[node->name];
        for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
        {
            if(sIter->second==specID)
            {
                (*cset)[sIter->second]=1;
            }
            else
            {
                (*cset)[sIter->second]=0;   
            }
        }
        int id=condsetMap_Tree.size();
        condsetMap_Tree[id]=cset;
        return 0;
    }
    //Otherwise get me all children under me and set their coordinates to 1
    map<string,int> leaves;
    getLeaves(node,leaves);
    cout <<"Internal node " << node->name << endl;
    for(map<string,int>::iterator nodeIter=leaves.begin();nodeIter!=leaves.end();nodeIter++)
    {
        cout <<nodeIter->first << endl;
    }
    map<int,int>* cset=new map<int,int>;
    for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
    {
        if(leaves.find(sIter->first)==leaves.end())
        {
            (*cset)[sIter->second]=0;
        }   
        else
        {
            (*cset)[sIter->second]=1;
        }
    }
    leaves.clear();
    int id=condsetMap_Tree.size();
    condsetMap_Tree[id]=cset;
    if(node->leftchild!=NULL)
    {
        initCondsetMap_Tree(node->leftchild);
    }
    if(node->rightchild!=NULL)
    {
        initCondsetMap_Tree(node->rightchild);
    }
    return 0;
}

int
MetaLearner::getLeaves(SpeciesDistance::Species* node, map<string,int>& leaves)
{
    if(node->leftchild==NULL && node->rightchild==NULL)
    {
        leaves[node->name]=0;
        return 0;
    }
    if(node->leftchild!=NULL)
    {
        getLeaves(node->leftchild,leaves);
    }
    if(node->rightchild!=NULL)
    {
        getLeaves(node->rightchild,leaves);
    }
    return 0;
}

int
MetaLearner::initCondsetMap_Nopool()
{
    int ind=0;
    for(map<int,string>::iterator eIter=speciesIDNameMap.begin();eIter!=speciesIDNameMap.end();eIter++)
    {
        INTINTMAP* cset=new INTINTMAP;

        for(map<int,string>::iterator dIter=speciesIDNameMap.begin();dIter!=speciesIDNameMap.end();dIter++)
        {
            if(eIter==dIter)
            {
                (*cset)[dIter->first]=1;
            }
            else
            {
                (*cset)[dIter->first]=0;
            }
        }
        int currind=(int)pow(2.0,ind);
        condsetMap[currind]=cset;
        string condKey;
        genCondSetKey(*cset,condKey);
        condsetKeyIDMap[condKey]=currind;
        condsetIDKeyMap[currind]=condKey;
        ind=ind+1;
    }
    return 0;
}

INTINTMAP*
MetaLearner::getConditionSet(int cind)
{
    if(condsetMap.find(cind)==condsetMap.end())
    {
        cout <<"Did not find any condition sets associated with " << cind << endl;
        exit(0);
    }
    return condsetMap[cind];
}


int
MetaLearner::collectMoves_Orthogroups(int currK)
{
    for(int i=0;i<moveSet.size();i++)
    {
        delete moveSet[i];
    }
    moveSet.clear();
    map<int,MappedOrthogroup*>& orthogroupSet=ogr->getMappedOrthogroups();

    //cout << "collectMoves_Orthogroups" << endl;
    cout << "inputOGList.size() = " << inputOGList.size() << endl;
    cout << "inputRegulatorOGs.size() = " << inputRegulatorOGs.size() << endl;

    //Now we will have a move for one orthogroup at a time
    for(map<int,MappedOrthogroup*>::iterator oIter=orthogroupSet.begin();oIter!=orthogroupSet.end();oIter++)
    {
        // cout << "  oIter->first = " << oIter->first << endl;
        if(inputOGList.size()>0  && inputOGList.find(oIter->first)==inputOGList.end())
        {
            continue;
        }
        /*if(oIter->first==189)
        {   
    //      cout <<"Stop here" << endl;
        }*/
        MappedOrthogroup* targetogrp=oIter->second;
        map<string,GeneMap*>& targetgrpMembers=targetogrp->getOrthoMembers();
        // cout << "  targetgrpMembers.size() = " << targetgrpMembers.size() << endl;
        map<int,double> bestscore_PerSpecies;
        map<int,double> bestscoreImprovement_PerSpecies;
        map<int,int> besttarget_PerSpecies;
        map<int,int> besttf_PerSpecies;
        //map<int,INTDBLMAP*> bestregWt_PerSpecies;
        int bestcsetid=-1;
        //The logic of this is we will basically search for the utility of each regulator across every species. The datalikelihood
        //term is computed separately from the prior. Then we will consider what will happen if were to make moves for all species.
        double bestScoreImprovement_TF=-1;
        double bestprior=0;
        for(map<int,int>::iterator regOGIter=inputRegulatorOGs.begin();regOGIter!=inputRegulatorOGs.end();regOGIter++)
        {
            // cout << "  regOGIter->first = " << regOGIter->first << endl;
            if(regOGIter->first==oIter->first)
            {
                continue;
            }
            /*if(regOGIter->first==26)
            {
                //cout << "Stop here" << endl;
            }*/
            map<int,double> score_PerSpecies;
            map<int,double> scoreImprovement_PerSpecies;
            map<int,int> target_PerSpecies;
            map<int,int> tf_PerSpecies;
            //map<int,INTDBLMAP*> regWt_PerSpecies;
            if(orthogroupSet.find(regOGIter->first)==orthogroupSet.end())
            {
                //cout << "Regulator ID " << regOGIter->first << " not found" << endl; 
                continue;
            }
            MappedOrthogroup* tfogrp=orthogroupSet[regOGIter->first];
            map<string,GeneMap*>& tfgrpMembers=tfogrp->getOrthoMembers();
            char ogPairKey[256];
            sprintf(ogPairKey,"%d-%d",regOGIter->first,oIter->first);
            string ogPair(ogPairKey);
            double oldpriorScore=ogpairPrior[ogPair];
            //for(map<string,GeneMap*>::iterator specIter=targetgrpMembers.begin();specIter!=targetgrpMembers.end();specIter++)
            for(map<string,SpeciesDataManager*>::iterator specIter=speciesDataSet.begin();specIter!=speciesDataSet.end();specIter++)
            {
                //Let's see if the regulator OG exists in this species
                if(tfgrpMembers.find(specIter->first)==tfgrpMembers.end())
                {
                    continue;
                }
                if(targetgrpMembers.find(specIter->first)==targetgrpMembers.end())
                {
                    continue;
                }
                /*if(speciesDataSet.find(specIter->first)==speciesDataSet.end())
                {
    //              continue;
                }*/
                int specID=speciesNameIDMap[specIter->first];
                SpeciesDataManager* sdm=speciesDataSet[specIter->first];
                map<string,int>& candidateregulators_Species=sdm->getRegulators();
                FactorGraph* speciesGraph=sdm->getFactorGraph();
                VariableManager* vMgr=sdm->getVariableManager();
                VSET& varSet=vMgr->getVariableSet();
                //GeneMap* geneMap_Tgts=specIter->second;
                GeneMap* geneMap_Tgts=targetgrpMembers[specIter->first];
                map<string,map<string,STRINTMAP*>*>& speciesTargetSet=geneMap_Tgts->getGeneSet();
                GeneMap* geneMap_TFs=tfgrpMembers[specIter->first];
                map<string,map<string,STRINTMAP*>*>& speciesTFSet=geneMap_TFs->getGeneSet();
                //If the species has multiple genes in it's list, we consider the member which will give the max improvement, that is
                //highest data likelihood. Similarly, we will extract the highest likelihood if there are multiple regulators
                double maxScore=-999999;
                double maxScoreImprovement=-999999;
                //Best here makes sense only if there are multiple TFs and targets
                int bestTarget=-1;
                int bestTF=-1;
                //INTDBLMAP wts;
                for(map<string,map<string,STRINTMAP*>*>::iterator vIter=speciesTargetSet.begin();vIter!=speciesTargetSet.end();vIter++)
                {   
                    int targetID=vMgr->getVarID(vIter->first.c_str());
                    if(targetID==-1)
                    {
                        continue;
                    }
                    Variable* target=varSet[targetID];
                    SlimFactor* sFactor=speciesGraph->getFactorAt(targetID);
                    for(map<string,map<string,STRINTMAP*>*>::iterator rIter=speciesTFSet.begin();rIter!=speciesTFSet.end();rIter++)
                    {
                        int regID=vMgr->getVarID(rIter->first.c_str());
                        if(regID==-1)
                        {
                            continue;
                        }
                        if(candidateregulators_Species.find(rIter->first)==candidateregulators_Species.end())
                        {
                            continue;
                        }
                        if(regID==targetID)
                        {
                            continue;
                        }
                        //If the edge already exists in the MB of sFactor continue  
                        if(sFactor->mergedMB.find(regID)!=sFactor->mergedMB.end())
                        {
                            continue;
                        }
                        if(!checkMBSize(specID,targetID,currK))
                        {
                            continue;
                        }
                        Variable* reg=varSet[regID];

                        /*if((strcmp(reg->getName().c_str(),"YBR289W")==0) && (strcmp(target->getName().c_str(),"YER037W")==0))
                        {
                //          cout <<"Stop here " << endl;
                        }*/
                        //Otherwise get the score delta of adding this reguluator in sFactor's MB.
                        double scoreImprovement=0;
                        double newScore=0;
                        double newTargetScore=0;
                        //INTDBLMAP regwt;
                        INTINTMAP* cset=getConditionSet(specID);
                        // shilu: no need to compute regwt at this time
                        //getNewPLLScore(specID,*cset,reg,target,newScore,scoreImprovement,oIter->first,regwt);
                        getNewPLLScore(specID,*cset,reg,target,newScore,scoreImprovement,oIter->first);
                        //cout << "specID=" << specID <<  " regOGID=" << regOGIter->first << " targetOGID=" << oIter->first << " TFname=" << rIter->first <<" targetgene=" << vIter->first << " regvaribleID=" << regID << " targetvaribleID=" <<targetID << " scoreImprovement=" << scoreImprovement << " newScore=" << newScore << endl;
                        //cout <<specIter->first << " " <<vIter->first << "-" <<rIter->first <<" scoreImprovement=" << scoreImprovement << endl;
                        if((scoreImprovement>0) && (scoreImprovement>maxScoreImprovement))
                        {
                            bestTarget=targetID;
                            bestTF=regID;
                            maxScoreImprovement=scoreImprovement;
                            maxScore=newScore;
                            /*wts.clear();
                            for(INTDBLMAP_ITER wIter=regwt.begin();wIter!=regwt.end();wIter++)
                            {
                                wts[wIter->first]=wIter->second;
                            }*/
                        }
                        //regwt.clear();
                    }
                }
                if(maxScoreImprovement<0)
                {
                    //wts.clear();
                    continue;
                }
                //At this stage we are done with this species, and if maxDLL >0 we proceed with updating the information for this species
                score_PerSpecies[specID]=maxScore;
                scoreImprovement_PerSpecies[specID]=maxScoreImprovement;
                target_PerSpecies[specID]=bestTarget;
                tf_PerSpecies[specID]=bestTF;
                /*INTDBLMAP* regwtforspecies=new INTDBLMAP;
                for(INTDBLMAP_ITER wIter=wts.begin();wIter!=wts.end();wIter++)
                {
                    (*regwtforspecies)[wIter->first]=wIter->second;
                }
                regWt_PerSpecies[specID]=regwtforspecies;*/
            }
            //At this stage we have considered the candidate regulator in the orthogroup regOGIter in all species for our input gene
            if(scoreImprovement_PerSpecies.size()==0)
            {
                continue;
                //cout <<"Skipping no hit regulator " <<  regOGIter->first << endl;
            }
            else
            {
                /*cout <<"Edge could be interesting in " ;
                for(map<int,double>::iterator sIter=scoreImprovement_PerSpecies.begin();sIter!=scoreImprovement_PerSpecies.end();sIter++)
                {
                    cout <<" " << speciesIDNameMap[sIter->first];
                }
                cout << endl;*/
            }
            double bestImprovement=0;
            int csetid=-1;
            double bestPrior_Reg=0;
            //cout <<"Compute for each condsetMap_Tree" << endl;
            //Now we wish to see how good it would be add these edges in different groups
            map<string,int> speciesEdgeStat;
            for(map<int,INTINTMAP*>::iterator setIter=condsetMap_Tree.begin();setIter!=condsetMap_Tree.end();setIter++)
            {
                INTINTMAP* cset=setIter->second;
                int specmatch=0;
                int expmatch=0;
                double netImprovement=0;
                speciesEdgeStat.clear();
                for(INTINTMAP_ITER cIter=cset->begin();cIter!=cset->end();cIter++)
                {
                    if(cIter->second==1)
                    {
                        expmatch++;
                        if(scoreImprovement_PerSpecies.find(cIter->first)==scoreImprovement_PerSpecies.end())
                        {
                            continue;
                        }
                        netImprovement=netImprovement+scoreImprovement_PerSpecies[cIter->first];
                        specmatch++;
                        speciesEdgeStat[speciesIDNameMap[cIter->first]]=1;
                    }
                    else
                    {
                        speciesEdgeStat[speciesIDNameMap[cIter->first]]=0;
                    }
                }
                /*if(scoreImprovement_PerSpecies.size()>1)
                {
                    if(specmatch==expmatch)
                    {
                    //  cout <<"Stop here " << endl;
                    }
                }*/
                if(specmatch<expmatch)
                {
                    //This configuration is not valid
                    continue;
                }
                double ePrior=speciesData->getEdgeStatusProb(speciesEdgeStat);
                ePrior=log(ePrior);
                netImprovement=netImprovement+ePrior-oldpriorScore;
                //cout << "condition"<<setIter->first <<": ePrior=" << ePrior << " oldpriorScore=" << oldpriorScore << " netImprovementwPrior=" << netImprovement << endl;
                if(netImprovement>bestImprovement)
                {
                    bestImprovement=netImprovement;
                    csetid=setIter->first;
                    bestPrior_Reg=ePrior;
                }
            }
            if(csetid==-1)
            {
                continue;
            }
            if(bestImprovement>bestScoreImprovement_TF)
            {
                INTINTMAP* cset=condsetMap_Tree[csetid];
                /*cout <<"OG" << oIter->first<<  " Updating  with regulator " << regOGIter->first << " with  species ";
                for(map<int,int>::iterator sIter=cset->begin();sIter!=cset->end();sIter++)
                {
                    if(sIter->second==0)
                    {
                        continue;
                    }
                    cout <<" " << speciesIDNameMap[sIter->first];
                }
                cout << endl;*/
                bestscore_PerSpecies.clear();
                bestscoreImprovement_PerSpecies.clear();
                besttarget_PerSpecies.clear();
                besttf_PerSpecies.clear();
                /*for(map<int,INTDBLMAP*>::iterator wtIter=bestregWt_PerSpecies.begin();wtIter!=bestregWt_PerSpecies.end();wtIter++)
                {
                    wtIter->second->clear();
                    delete wtIter->second;
                }
                bestregWt_PerSpecies.clear();*/
                for(INTINTMAP_ITER csIter=cset->begin();csIter!=cset->end();csIter++)
                {
                    if(csIter->second==0)
                    {
                        continue;
                    }
                    bestscore_PerSpecies[csIter->first]=score_PerSpecies[csIter->first];
                    bestscoreImprovement_PerSpecies[csIter->first]=scoreImprovement_PerSpecies[csIter->first];
                    besttarget_PerSpecies[csIter->first]=target_PerSpecies[csIter->first];
                    besttf_PerSpecies[csIter->first]=tf_PerSpecies[csIter->first];
                    //bestregWt_PerSpecies[csIter->first]=regWt_PerSpecies[csIter->first];
                }
                bestcsetid=csetid;
                bestprior=bestPrior_Reg;
                bestScoreImprovement_TF=bestImprovement;
            }
            /*else
            {
                for(map<int,INTDBLMAP*>::iterator wtIter=regWt_PerSpecies.begin();wtIter!=regWt_PerSpecies.end();wtIter++)
                {
                    wtIter->second->clear();
                    delete wtIter->second;
                }
                regWt_PerSpecies.clear();
            }*/
        }
        if(bestcsetid<0)
        {
            continue;
        }
        INTINTMAP* cset=condsetMap_Tree[bestcsetid];
        for(INTINTMAP_ITER csIter=cset->begin();csIter!=cset->end();csIter++)
        {
            if(csIter->second==0)
            {
                continue;
            }   
            MetaMove* move=new MetaMove;
            int tfid=besttf_PerSpecies[csIter->first];
            int tgtid=besttarget_PerSpecies[csIter->first];
            move->setSrcVertex(tfid);
            move->setConditionSetInd(csIter->first);
            SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[csIter->first]];
            VariableManager* vMgr=sdm->getVariableManager();
            VSET& varSet=vMgr->getVariableSet();
            //cout <<"Found edge for "<< speciesIDNameMap[csIter->first] <<" TFvar=" << tfid << " Targetvar=" << tgtid <<  " reg=" << varSet[tfid]->getName().c_str()  << " targeti="<< varSet[tgtid]->getName().c_str() << " score deltas=" << bestscoreImprovement_PerSpecies[csIter->first] <<endl;
            move->setTargetVertex(besttarget_PerSpecies[csIter->first]);
            move->setTargetMBScore(bestscore_PerSpecies[csIter->first]);
            move->setScoreImprovement(bestscoreImprovement_PerSpecies[csIter->first]);
            /*INTDBLMAP* wts=bestregWt_PerSpecies[csIter->first];
            if(wts->find(tfid)==wts->end())
            {
                cout <<"Did not find " << varSet[tfid]->getName().c_str() << " regulator in " << varSet[tgtid]->getName() << endl;
            }
            move->setSrcWeight(*bestregWt_PerSpecies[csIter->first]);*/
            moveSet.push_back(move);
        }
        besttarget_PerSpecies.clear();
        bestscore_PerSpecies.clear();
        //bestregWt_PerSpecies.clear();
        besttf_PerSpecies.clear();
        bestscoreImprovement_PerSpecies.clear();
    }

    return 0;
}


/*
int
MetaLearner::collectMoves_Deletions()
{
    return 0;
    for(int i=0;i<moveSet.size();i++)
    {
        delete moveSet[i];
    }
    moveSet.clear();
    //Go over all the edges right now and figure out if we want to remove any
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        SpeciesDataManager* sdm=sIter->second;
        FactorGraph* factorGraph=sdm->getFactorGraph();
        VariableManager* varMgr=sdm->getVariableManager();
        VSET& varSet=varMgr->getVariableSet();
        map<int,SlimFactor*>& factorSet=factorGraph->getAllFactors();
        for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
        {
            SlimFactor* sFactor=fIter->second;
            Variable* target=varSet[fIter->first];
            INTINTMAP& mbVars=sFactor->mergedMB;
            INTINTMAP saveMBVars;
            for(INTINTMAP_ITER mIter=mbVars.begin();mIter!=mbVars.end();mIter++)
            {
                saveMBVars[mIter->first]=mIter->second;
            }
            double maxScoreImpr=0;
            int worseRegulatorID=0;
            double bestScore=0;
            INTDBLMAP worstRegWts;
            if(saveMBVars.size()>1)
            {
                cout <<"Stop here" << endl;
            }
            for(INTINTMAP_ITER mIter=saveMBVars.begin();mIter!=saveMBVars.end();mIter++)
            {
                Variable* regulator=varSet[mIter->first];
                INTINTMAP_ITER vIter=sFactor->mergedMB.find(mIter->first);
                sFactor->mergedMB.erase(vIter);
                double scoreImprovement=0;
                double newScore=0;
                double newTargetScore=0;
                int status=0;
                double priorScore=0;
                double oldpriorScore=0;
                if(speciesData->getRoot()!=NULL)
                {
                    int uogno=ogr->getMappedOrthogroupID(regulator->getName().c_str(),sIter->first.c_str());
                    int vogno=ogr->getMappedOrthogroupID(target->getName().c_str(),sIter->first.c_str());
                    if(vogno!=-1)
                    {   
                        priorScore=getEdgePrior(speciesNameIDMap[sIter->first],regulator,target,vogno,0);
                        char ogPairKey[256];
                        sprintf(ogPairKey,"%d-%d",uogno,vogno);
                        string ogPair(ogPairKey);
                        oldpriorScore=ogpairPrior[ogPair];
                    }
                }
                INTDBLMAP regwt;
                double pll_d=getPLLScore_Condition_Tracetrick((string&)sIter->first,sFactor,status,regwt);
                //double pll_d=getPLLScore_Condition((string&)sIter->first,sFactor,status,regwt);
                //double pll_d=getPLLScore_Condition((string&)sIter->first,sFactor);
                double dImpr=pll_d+priorScore-oldpriorScore;
                if(dImpr>maxScoreImpr)
                {
                    worseRegulatorID=mIter->first;
                    maxScoreImpr=dImpr;
                    worstRegWts.clear();
                    for(INTDBLMAP_ITER wIter=regwt.begin();wIter!=regwt.end();wIter++)
                    {
                        worstRegWts[wIter->first]=wIter->second;
                    }
                }
                regwt.clear();
                sFactor->mergedMB[mIter->first]=0;
            }
            if(maxScoreImpr<=0)
            {
                worstRegWts.clear();
                continue;
            }
            MetaMove* move=new MetaMove;
            move->setSrcVertex(worseRegulatorID);
            move->setConditionSetInd(speciesNameIDMap[sIter->first]);
            cout <<" Found edge for deletion " <<varSet[worseRegulatorID]->getName().c_str() << "<->" 
            << varSet[fIter->first]->getName() << " score deltas" << maxScoreImpr << endl;
            move->setTargetVertex(fIter->first);
            move->setTargetMBScore(bestScore);
            move->setScoreImprovement(maxScoreImpr);
            move->setSrcWeight(worstRegWts);
            moveSet.push_back(move);
        }
    }
    return 0;
}*/


int
MetaLearner::getNewPLLScore(int cid, INTINTMAP& conditionSet, Variable* u, Variable* v, double& targetmbScore, double& scoreImprovement,int orthoGrpNo) //INTDBLMAP& regwt
{
    //Each condition set has an evidence manager which pools the data from the specific conditions.
    string condKey;
    genCondSetKey(conditionSet,condKey);
    SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[cid]];
    PotentialManager* potMgr=sdm->getPotentialManager();
    VSET& varSet=sdm->getVariableManager()->getVariableSet();
    FactorGraph* fg=sdm->getFactorGraph();
    SlimFactor* sFactor=fg->getFactorAt(u->getID());
    SlimFactor* dFactor=fg->getFactorAt(v->getID());
    /*if(strcmp(u->getName().c_str(),"YMR037C")==0)
    {
        cout <<"Stop here "<< endl;
    }*/
    map<int,double>* varNeighborhoodPrior=varNeighborhoodPrior_PerSpecies[speciesIDNameMap[cid]];
    map<string,double>* edgePresenceProb=edgePresenceProb_PerSpecies[speciesIDNameMap[cid]];
    double currPrior=(*varNeighborhoodPrior)[v->getID()];
    bool toDel_d=true;
    if(dFactor->mergedMB.find(u->getID())!=dFactor->mergedMB.end())
    {
        toDel_d=false;
    }
    double plus=0;
    double minus=0;
    dFactor->mergedMB.insert(u->getID()); //Aug 23: dFactor->mergedMB[u->getID()]=0;
    int status=0;
    for(auto mIter=dFactor->mergedMB.begin();mIter!=dFactor->mergedMB.end();mIter++)
    {
        Variable* aVar=varSet[*mIter];
        string regulatorKey(aVar->getName().c_str());
        regulatorKey.append("\t");
        regulatorKey.append(v->getName().c_str());
        double p=(*edgePresenceProb)[regulatorKey];
        minus=minus+log(1-p);
        plus=plus+log(p);
    }
    double pll_d=getPLLScore_Condition_Tracetrick(speciesIDNameMap[cid],dFactor,status);
    //double pll_d=getPLLScore_Condition_Tracetrick(speciesIDNameMap[cid],dFactor,status,regwt);
    //double pll_d=getPLLScore_Condition(speciesIDNameMap[cid],dFactor,status,regwt);
    if(status==-1)
    {
        scoreImprovement=-1;
        if(toDel_d)
        {
            auto dIter=dFactor->mergedMB.find(u->getID());
            dFactor->mergedMB.erase(dIter);
        }
        return 0;
    }
    currPrior=currPrior+plus-minus;
    pll_d=pll_d+currPrior;
    /*double priorScore=0;
    double oldpriorScore=0;
    if(speciesData->getRoot()!=NULL)
    {
        priorScore=getEdgePrior(cid,u,v,orthoGrpNo);
        int uogno=ogr->getMappedOrthogroupID(u->getName().c_str(),speciesIDNameMap[cid].c_str());
        int vogno=ogr->getMappedOrthogroupID(v->getName().c_str(),speciesIDNameMap[cid].c_str());
        if(vogno!=-1)
        {
            char ogPairKey[256];
            sprintf(ogPairKey,"%d-%d",uogno,vogno);
            string ogPair(ogPairKey);
            oldpriorScore=ogpairPrior[ogPair];
        }
    }*/
    targetmbScore=pll_d;
    scoreImprovement=0;
    //double dImpr=(targetmbScore+priorScore)-(dFactor->mbScore+oldpriorScore);
    //Don't include the prior. Just use the data likelihood improvement
    double dImpr=targetmbScore-dFactor->mbScore;
    if(dImpr<=0)
    {
        scoreImprovement=-1;
    }
    else
    {
        scoreImprovement=dImpr;
    }
    if(toDel_d)
    {
        auto dIter=dFactor->mergedMB.find(u->getID());
        dFactor->mergedMB.erase(dIter);
    }
    
    return 0;
}


double
MetaLearner::getPLLScore_Condition(string& speciesName,SlimFactor* sFactor)
{
    SpeciesDataManager* spd=speciesDataSet[speciesName];
    PotentialManager* potMgr=spd->getPotentialManager();
    EvidenceManager* evMgr=spd->getEvidenceManager();
    VariableManager* varMgr=spd->getVariableManager();
    VSET& varSet=varMgr->getVariableSet();  
    /*string mbkey;
    char mbchar[256];
    for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
    {
        sprintf(mbchar,"-%d",mIter->first);
        mbkey.append(mbchar);
    }
    
    if(strcmp(mbkey.c_str(),"-1-35-39-45")==0)
    {
        cout <<"Found " << mbkey.c_str() << endl;
    }*/
    double unreg_pll=potMgr->getPseudoLikelihood(sFactor,varSet,false);
    if((isnan(unreg_pll)) || (isinf(unreg_pll)))
    {
        cout <<"Found nan for " << sFactor->fId<< ": MB: ";
        for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
        {
            cout <<"-"<<*mIter;
        }
        cout<< endl;
    }
    /* //comment by shilu
    //double varCnt=sFactor->mergedMB.size()+1;
    double varCnt=sFactor->mergedMB.size();
    double paramCnt=2*varCnt;
    //double paramCnt=varCnt;
    paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
    double pll=unreg_pll-((paramCnt/2)*log(evMgr->getTestSet().size()));
    //pll=unreg_pll-(0.5*log(evMgr->getTestSet().size()));
    pll=unreg_pll;*/
    return unreg_pll;
}

/*
// commented by vperiyasamy
double
MetaLearner::getPLLScore_Condition(string& speciesName,SlimFactor* sFactor,int& status,INTDBLMAP& regWts)
{
    PotentialManager* potMgr=speciesDataSet[speciesName]->getPotentialManager(); // grab managers
    EvidenceManager* evMgr=speciesDataSet[speciesName]->getEvidenceManager(); // for this
    VariableManager* varMgr=speciesDataSet[speciesName]->getVariableManager(); // species
    VSET& varSet=varMgr->getVariableSet();  

    string mbkey;
    char mbchar[256];
    if(sFactor->mergedMB.size()>0) // sFactor here is for vertex v, if no parents then what?
    {
    //  cout <<"Stop here" << endl;
    }
    int statuspot=0;
    Potential* apot=new Potential;

    // get unregularized PLL score?
    double unreg_pll=potMgr->getPseudoLikelihood(sFactor,varSet,false,statuspot,apot);
    if(statuspot==-1)
    {
        status=-1;
        delete apot;
        return 0;
    }

    // get condition weights from filled in Potential
    INTDBLMAP& potwts=apot->getCondWeight();
    for(INTDBLMAP_ITER wIter=potwts.begin();wIter!=potwts.end();wIter++)
    {
        regWts[wIter->first]=wIter->second;
    }
        
    if((isnan(unreg_pll)) || (isinf(unreg_pll)))
    {
        cout <<"Found nan for " << sFactor->fId<< ": MB: ";
        for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
        {
            cout <<"-"<<mIter->first;
        }
        cout<< endl;
    }

    // increase parent count because we have a new edge 
    double varCnt=sFactor->mergedMB.size()+1.0;
    //double varCnt=sFactor->mergedMB.size();

    double paramCnt=2*varCnt;
    paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
    //paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
    double pll=unreg_pll-(0.5*paramCnt*log(evMgr->getTestSet().size())); // MERLIN uses training set size? does this make a difference
    //double pll=unreg_pll;
    pll=unreg_pll; // ?????? why are we just putting it back to unreg_pll
    delete apot;
    return pll;
}*/

double
//MetaLearner::getNewPLLScore_Condition_Tracetrick(int csetId, int vId, int uId, Potential* newPot)
MetaLearner::getPLLScore_Condition_Tracetrick(string& speciesName,SlimFactor* sFactor,int& status)
{
    PotentialManager* potMgr=speciesDataSet[speciesName]->getPotentialManager(); // grab managers
    EvidenceManager* evMgr=speciesDataSet[speciesName]->getEvidenceManager(); // for this
    //VariableManager* varMgr=speciesDataSet[speciesName]->getVariableManager(); // species
    //VSET& varSet=varMgr->getVariableSet();
    //cout << speciesName <<  " sFactor=" <<sFactor->fId << " "<< endl;
    //vector<int> * tSet = &evMgr->getTestSet();
    //int datasize = tSet->size();
    int datasize = evMgr->getTestSetSize(); 
    double pll=potMgr->computePotentialMBCovMean(sFactor,status,datasize);
    //Need to fix this to be set automatically
    //get parameter prior
    //double paramPrior=0;
    //cout << speciesName <<  " sFactor=" <<sFactor->fId << " pll=" << pll << endl;
    return pll;
}

/*
// added by vperiyasamy
double
//MetaLearner::getNewPLLScore_Condition_Tracetrick(int csetId, int vId, int uId, Potential* newPot)
MetaLearner::getPLLScore_Condition_Tracetrick(string& speciesName,SlimFactor* sFactor,int& status,INTDBLMAP& regWts)
{
    PotentialManager* potMgr=speciesDataSet[speciesName]->getPotentialManager(); // grab managers
    EvidenceManager* evMgr=speciesDataSet[speciesName]->getEvidenceManager(); // for this
    VariableManager* varMgr=speciesDataSet[speciesName]->getVariableManager(); // species
    VSET& varSet=varMgr->getVariableSet();  

    string mbkey;
    char mbchar[256];
    if(sFactor->mergedMB.size()>0) { // sFactor here is for vertex v, if no parents then what?
    //  cout <<"Stop here" << endl;
    }

    // first define newPot and set it up
    Potential* newPot = new Potential;
    newPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR); // set v as the target for this potential

    for(INTINTMAP_ITER aIter=sFactor->mergedMB.begin();aIter!=sFactor->mergedMB.end();aIter++) {
        Variable* aVar=varSet[aIter->first];
        newPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
    }
    newPot->potZeroInit();

    // DO WE NEED THIS?
    //newPot->setCondBias(sFactor->potFunc->getCondBias());
    //newPot->setCondVariance(sFactor->potFunc->getCondVariance());
    //newPot->setCondWeight(sFactor->potFunc->getCondWeight());

    potMgr->populatePotential(newPot,false);
    //This function creates a submatrix of the covariance matrix and inverts it
    newPot->initMBCovMean();
    if(newPot->getCondVariance() < 0) {
        status=-1;
        delete newPot;
        return 0;
    }

    // get condition weights from filled in Potential
    INTDBLMAP& potwts = newPot->getCondWeight();
    for(INTDBLMAP_ITER wIter = potwts.begin(); wIter != potwts.end(); wIter++) {
        regWts[wIter->first] = wIter->second;
    }

    double pll=0; 
    //Need to fix this to be set automatically
    //get parameter prior
    double paramPrior=0;

    // now define parentPot and set it up

    Potential* parentPot=new Potential;
    VSET& vars=newPot->getAssocVariables();
    for(VSET_ITER vIter=vars.begin();vIter!=vars.end();vIter++)
    {
        if(vIter->first==sFactor->fId)
        {
            continue;
        }
        Variable* aVar=vars[vIter->first];
        parentPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
    }
    parentPot->potZeroInit();
    potMgr->populatePotential(parentPot,false);
    //parentPot->initMBCovMean();

    // NOTE: use test set in order to get same results as non tracetrick
    INTINTMAP* tSet = &evMgr->getTestSet();
    int datasize = tSet->size();
    double jointll1 = newPot->computeLL_Tracetrick(datasize);
    double jointll2 = parentPot->computeLL_Tracetrick(datasize);
    //double vCnt=(double)newPot->getAssocVariables().size();
    //double paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
    pll = jointll1 - jointll2;
    //pll=pll-(lambda*paramCnt*log(datasize));
    //pll = pll - (0.5 * paramCnt * log(datasize));
    delete newPot;
    delete parentPot;
    return pll;
}
*/
double
MetaLearner::getValidationPLLScore_Condition(int cId, int vId)
{
    SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[cId]];
    EvidenceManager* evMgr=sdm->getEvidenceManager();
    VariableManager* varMgr=sdm->getVariableManager();
    INTINTMAP& vSet=evMgr->getValidationSet();
    double pll=0; 
    //Need to fix this to be set automatically
    double wt=0.5;
    if(!dataPool)
    {
        wt=1;
    }
    double paramCnt=0;
    VSET& varSet=varMgr->getVariableSet();
    map<int,Potential*> potSet;
    int thresholded=0;
    for(INTINTMAP_ITER eIter=vSet.begin();eIter!=vSet.end();eIter++)
    {
        EMAP* evidMap=evMgr->getEvidenceAt(eIter->first);
        //Go over all condition sets that include cInd
        double cll=0;
        for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
        {
            INTINTMAP* cset=csIter->second;
            if((*cset)[cId]==0)
            {
                continue;
            }
            SpeciesDataManager* localSdm=speciesDataSet[speciesIDNameMap[csIter->first]];
            FactorGraph* fg=localSdm->getFactorGraph();
            SlimFactor* sFactor=fg->getFactorAt(vId);
            Potential* sPot=NULL;
            if(potSet.find(csIter->first)==potSet.end())
            {
                sPot=new Potential;
                potSet[csIter->first]=sPot;
                sPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
                for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
                {
                    Variable* aVar=varSet[*mIter];
                    sPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
                }
                sPot->potZeroInit();
                string condKey;
                genCondSetKey(*cset,condKey);
                PotentialManager* potMgr=pooledPotentials[condKey];
                //potMgr->populatePotential(sPot,false);
                potMgr->populatePotential_Eff(sPot,false);
                sPot->initMBCovMean();
            }
            else
            {
                sPot=potSet[csIter->first];
            }
            double pval=sPot->getCondPotValueFor(evidMap);

            if(pval<1e-50)
            {
                pval=1e-50;
                thresholded++;
            }
            if((pval==0) || (isnan(pval)) || (isinf(pval)))
            {
            //  cout <<"Stop here" << endl;
            }
            cll=cll+(pval*wt);
            if(eIter==vSet.begin())
            {
                double vCnt=(double)sPot->getAssocVariables().size();
                paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
            }
        }
        //Check here for really small loglikelihoods
        pll=pll+log(cll);

    }
    pll=pll-(0.5*paramCnt*log(vSet.size()));
    for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
    {
        delete pIter->second;
    }
    potSet.clear();
    return pll;
}

//This function computes the prior score of the edge between u and v in the species specified by conditionset
double
MetaLearner::getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo)
{
    double prior=0;
    if(speciesDataSet.size()==1)
    {
        return 0;
    }
    STRINTMAP edgeStatus;
    string& species=speciesIDNameMap[specID];
    for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
    {
        if(sIter->second==specID)
        {
            continue;
        }
        int hit=0;
        STRINTMAP* specTFs=ogr->getOrtholog(species.c_str(),u->getName().c_str(),sIter->first.c_str());
        STRINTMAP* specTargets=ogr->getOrtholog(species.c_str(),v->getName().c_str(),sIter->first.c_str());
        if((specTFs==NULL) || (specTargets==NULL))
        {
            edgeStatus[sIter->first]=0;
        }   
        else
        {
            SpeciesDataManager* sdm=speciesDataSet[sIter->first];
            FactorGraph* speciesGraph=sdm->getFactorGraph();
            VariableManager* vMgr=sdm->getVariableManager();
            VSET& varSet=vMgr->getVariableSet();
            //Do any of the genes in this species have data and if is anybody connected to u already
            for(STRINTMAP_ITER gIter=specTargets->begin();gIter!=specTargets->end();gIter++)
            {
                int tgtID=vMgr->getVarID(gIter->first.c_str());
                if(tgtID==-1)
                {
                    continue;
                }
                SlimFactor* sFactor=speciesGraph->getFactorAt(tgtID);
                for(STRINTMAP_ITER fIter=specTFs->begin();fIter!=specTFs->end();fIter++)
                {
                    int tfID=vMgr->getVarID(fIter->first.c_str());
                    if(tfID==-1)
                    {
                        continue;
                    }
                    if(sFactor->mergedMB.find(tfID)!=sFactor->mergedMB.end())
                    {
                        hit++;
                    }
                }
            }
        }
        if(hit>0)
        {
            edgeStatus[sIter->first]=1;
        }
        else
        {
            edgeStatus[sIter->first]=0;;
        }
    }
    edgeStatus[species]=1;
    prior=speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior=log(prior);
    edgeStatus.clear();
    return logPrior;
}


double 
MetaLearner::getEdgePrior_PerSpecies(int tfID, int targetID,SpeciesDataManager* sdm)
{
    INTDBLMAP* regPriors=NULL;
    double prior=1/(1+exp(-1*beta1));
    double motifweight=0;
    double chipweight=0;
    map<int,map<int,double>*>& motifNetwork=sdm->getMotifNetwork();
    if(motifNetwork.find(tfID)!=motifNetwork.end())
    {
        map<int,double>* values=motifNetwork[tfID];
        if(values->find(targetID)!=values->end())
        {
            motifweight=(*values)[targetID];
        }
    }/*else{
        cout << "tfID=" <<tfID << " targetID=" << targetID << " not found" <<endl;
    }*/
    /*if(edgePriors_Motif.find(targetID)!=edgePriors_Motif.end())
    {
        regPriors=edgePriors_Motif[targetID];
        if(regPriors->find(tfID)!=regPriors->end())
        {
            motifweight=(*regPriors)[tfID];
        }
    }
    if(edgePriors_ChIP.find(targetID)!=edgePriors_ChIP.end())
    {   
        regPriors=edgePriors_ChIP[targetID];
        if(regPriors->find(tfID)!=regPriors->end())
        {
            chipweight=(*regPriors)[tfID];
        }
    }
    double fwt=(chipweight*beta_chip) + (motifweight*beta_chip);*/
    double fwt=motifweight*beta2;
    prior=1/(1+exp(-1*(beta1+fwt)));
    //cout << "motifweight=" << motifweight<<" fwt=" << fwt << " prior=" << prior << endl;
    if(prior<1e-6)
    {
        prior=1e-6;
    }
    if(prior==1)
    {
        prior=1-1e-6;
    }
    return prior;
}

//This function computes the prior score of the edge between u and v in the species specified by conditionset
double
MetaLearner::getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo, int addRemStat)
{
    double prior=0;
    STRINTMAP edgeStatus;
    string& species=speciesIDNameMap[specID];
    for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
    {
        int hit=0;
        STRINTMAP* specTargets=ogr->getOrtholog(species.c_str(),v->getName().c_str(),sIter->first.c_str());
        STRINTMAP* specTFs=ogr->getOrtholog(species.c_str(),u->getName().c_str(),sIter->first.c_str());
        if((specTargets==NULL)||(specTFs==NULL))
        {
            edgeStatus[sIter->first]=0;
        }   
        else
        {
            SpeciesDataManager* sdm=speciesDataSet[sIter->first];
            FactorGraph* speciesGraph=sdm->getFactorGraph();
            VariableManager* vMgr=sdm->getVariableManager();
            VSET& varSet=vMgr->getVariableSet();
            //Do any of the genes in this species have data and if is anybody connected to u already
            for(STRINTMAP_ITER gIter=specTargets->begin();gIter!=specTargets->end();gIter++)
            {
                int tgtID=vMgr->getVarID(gIter->first.c_str());
                if(tgtID==-1)
                {
                    continue;
                }
                SlimFactor* sFactor=speciesGraph->getFactorAt(tgtID);
                for(STRINTMAP_ITER fIter=specTFs->begin();fIter!=specTFs->end();fIter++)
                {
                    int tfID=vMgr->getVarID(fIter->first.c_str());
                    if(tfID==-1)
                    {
                        continue;
                    }
                    if(sFactor->mergedMB.find(tfID)!=sFactor->mergedMB.end())
                    {
                        hit++;
                    }
                }
            }
        }
        if(hit>0)
        {
            edgeStatus[sIter->first]=1;
        }
    }
    //Now we need to figure out what the status of the edge is if v has duplicates and it is present in the dataset. Also if u has
    //duplicates and is present in the dataset.
    SpeciesDataManager* sdm=speciesDataSet[species];
    FactorGraph* speciesGraph=sdm->getFactorGraph();
    VariableManager* vMgr=sdm->getVariableManager();
    MappedOrthogroup* ogrp_TF=ogr->getMappedOrthogroup(u->getName().c_str(),species.c_str());
    MappedOrthogroup* ogrp_Tgt=ogr->getMappedOrthogroup(v->getName().c_str(),species.c_str());
    map<string,map<string,STRINTMAP*>*>& targets=ogrp_Tgt->getSpeciesHits(species.c_str())->getGeneSet();
    map<string,map<string,STRINTMAP*>*>& tfs=ogrp_TF->getSpeciesHits(species.c_str())->getGeneSet();
    int hitinspecies=0;
    for(map<string,map<string,STRINTMAP*>*>::iterator tgtIter=targets.begin();tgtIter!=targets.end();tgtIter++)
    {
        int vId=vMgr->getVarID(tgtIter->first.c_str());
        if(vId==-1)
        {
            continue;
        }
        SlimFactor* sFactor=speciesGraph->getFactorAt(vId);
        for(map<string,map<string,STRINTMAP*>*>::iterator tfIter=tfs.begin();tfIter!=tfs.end();tfIter++)
        {
            int uId=vMgr->getVarID(tfIter->first.c_str());
            if(uId==-1)
            {
                continue;
            }
            if(sFactor->mergedMB.find(uId)!=sFactor->mergedMB.end())
            {
                hitinspecies++;
            }
        }
    }
    if(hitinspecies>0)
    {
        edgeStatus[species]=1;
    }
    else
    {
        edgeStatus[species]=0;
    }
    prior=speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior=log(prior);
    return logPrior;
}


int
MetaLearner::genCondSetKey(INTINTMAP& condSet, string& aKey)
{   
    char keypair[256];
    for(INTINTMAP_ITER cIter=condSet.begin();cIter!=condSet.end();cIter++)
    {
        sprintf(keypair,"-%d=%d",cIter->first,cIter->second);
        aKey.append(keypair);
    }
    return 0;
}

int 
MetaLearner::sortMoves()
{
    for(int m=0;m<moveSet.size();m++)
    {
        for(int n=m+1;n<moveSet.size();n++)
        {
            MetaMove* m1=moveSet[m];
            MetaMove* m2=moveSet[n];
            if(m1->getScoreImprovement()<m2->getScoreImprovement())
            {
                moveSet[m]=m2;
                moveSet[n]=m1;
            }
        }
    }
    return 0;
}

double
MetaLearner::makeMoves()
{
    map<int,INTINTMAP*> affectedVariables;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {
        int cind=speciesNameIDMap[gIter->first];
        INTINTMAP* csVars=new INTINTMAP;
        affectedVariables[cind]=csVars;
    }
    int successMove=0;
    double netScoreDelta=0;
    int net1Move=0;
    int net2Move=0;
    for(int m=0;m<moveSet.size();m++)
    {
        MetaMove* move=moveSet[m];
        if(attemptMove(move,affectedVariables)==0)
        {
            if(move->getConditionSetInd()==1)
            {
                net1Move++;
            }
            else if(move->getConditionSetInd()==2)
            {
                net2Move++;
            }
            successMove++;
            netScoreDelta=netScoreDelta+move->getScoreImprovement();
        }
        else
        {
            cout <<"Move not valid " << endl;
        }
    }
    for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
    {
        cIter->second->clear();
        delete cIter->second;
    }
    affectedVariables.clear();
    cout <<"Total successful moves " << successMove << " out of total " 
    << moveSet.size() << " with net score improvement " 
    << netScoreDelta
    << " net1 moves " << net1Move
    << " net2 moves " << net2Move
    << endl;
    return netScoreDelta;
}


/*int
MetaLearner::makeDelMoves()
{
    map<int,INTINTMAP*> affectedVariables;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {
        int cind=speciesNameIDMap[gIter->first];
        INTINTMAP* csVars=new INTINTMAP;
        affectedVariables[cind]=csVars;
    }
    int successMove=0;
    double netScoreDelta=0;
    int net1Move=0;
    int net2Move=0;
    for(int m=0;m<moveSet.size();m++)
    {
        MetaMove* move=moveSet[m];
        if(attemptDelMove(move,affectedVariables)==0)
        {
            successMove++;
            netScoreDelta=netScoreDelta+move->getScoreImprovement();
        }
    }
    for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
    {
        cIter->second->clear();
        delete cIter->second;
    }
    affectedVariables.clear();
    cout <<"Total successful moves " << successMove << " out of total " 
    << moveSet.size() << " with net score improvement " 
    << netScoreDelta
    << endl;
    return 0;
}*/

int
MetaLearner::attemptMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars)
{
    int specID=move->getConditionSetInd();
    SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[specID]];
    VSET& varSet=sdm->getVariableManager()->getVariableSet();
    Variable* u=varSet[move->getSrcVertex()];
    Variable* v=varSet[move->getTargetVertex()];
    int uogno=ogr->getMappedOrthogroupID(u->getName().c_str(),speciesIDNameMap[specID].c_str());
    int vogno=ogr->getMappedOrthogroupID(v->getName().c_str(),speciesIDNameMap[specID].c_str());

    INTINTMAP* csVars=affectedVars[specID];
    //INTINTMAP* csVars=affectedVars[1];

    /*if(csVars->find(uogno)!=csVars->end())
    {
        cout <<"TF OG " << uogno  << " already used for "<< speciesIDNameMap[specID] << endl;
        return -1;
    } */
    if(vogno==-1)
    {
        return 0;
    }
    if((vogno!=-1) && (csVars->find(vogno)!=csVars->end()))
    {
        cout <<"Target OG " << vogno  << " already used for "<< speciesIDNameMap[specID] << endl;
        return -1;
    }
    (*csVars)[vogno]=0;
    //(*csVars)[uogno]=0;
    FactorGraph* csGraph=sdm->getFactorGraph();
    SlimFactor* dFactor=csGraph->getFactorAt(move->getTargetVertex());
    dFactor->mergedMB.insert(move->getSrcVertex()); //dFactor->mergedMB[move->getSrcVertex()]=0;
    dFactor->mbScore=move->getTargetMBScore();
    //dFactor->setMBWts(move->getSrcWeight());
    char ogpair[256];
    sprintf(ogpair,"%d-%d",uogno,vogno);
    string ogpairKey(ogpair);
    //cout <<"Attempting move " << ogpairKey << " in species " << speciesIDNameMap[specID] << endl;
    if(edgeConditionMap.find(ogpairKey)==edgeConditionMap.end())
    {
        cout <<"No edge pair " << ogpairKey << endl;
        exit(0);
    }
    STRINTMAP* newEdgeStatus=NULL;
    if(affectedOGPairs.find(ogpairKey)==affectedOGPairs.end())
    {
        newEdgeStatus=new STRINTMAP;
        affectedOGPairs[ogpairKey]=newEdgeStatus;
    }
    else
    {
        newEdgeStatus=affectedOGPairs[ogpairKey];
    }
    string& species=speciesIDNameMap[specID];   
    //(*newEdgeStatus)[species]=1;
    (*newEdgeStatus)[species]=(*newEdgeStatus)[species]+1;
    if((*newEdgeStatus)[species]>1)
    {
        cout <<"Adding edge to same family again for "<< species << " gene " << v->getName() << " ogkey " << ogpairKey<< endl;
    }
    return 0;
}

/*
int
MetaLearner::attemptDelMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars)
{
    int specID=move->getConditionSetInd();
    SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[specID]];
    VSET& varSet=sdm->getVariableManager()->getVariableSet();
    Variable* u=varSet[move->getSrcVertex()];
    Variable* v=varSet[move->getTargetVertex()];
    int uogno=ogr->getMappedOrthogroupID(u->getName().c_str(),speciesIDNameMap[specID].c_str());
    int vogno=ogr->getMappedOrthogroupID(v->getName().c_str(),speciesIDNameMap[specID].c_str());
    
    INTINTMAP* csVars=affectedVars[1];
    if(csVars->find(uogno)!=csVars->end())
    {
        return -1;
    } 
    if((vogno!=-1) && (csVars->find(vogno)!=csVars->end()))
    {
        return -1;
    }
    FactorGraph* csGraph=sdm->getFactorGraph();
    SlimFactor* dFactor=csGraph->getFactorAt(move->getTargetVertex());
    INTINTMAP_ITER vIter=dFactor->mergedMB.find(move->getSrcVertex());
    if(vIter==dFactor->mergedMB.end())
    {
        cout <<"No mb var " << move->getSrcVertex() << " in mb of " << move->getTargetVertex() << endl;
        return 0;
    }
    (*csVars)[uogno]=0;
    dFactor->mergedMB.erase(vIter->first);
    dFactor->mbScore=move->getTargetMBScore();
    dFactor->setMBWts(move->getSrcWeight());
    if(vogno=-1)
    {
        return 0;
    }
    (*csVars)[vogno]=0;
    char ogpair[256];
    sprintf(ogpair,"%d-%d",uogno,vogno);
    string ogpairKey(ogpair);
    if(edgeConditionMap.find(ogpairKey)==edgeConditionMap.end())
    {
        cout <<"No edge pair " << ogpairKey << endl;
        exit(0);
    }
    INTINTMAP* currEdgeStatus=edgeConditionMap[ogpairKey];
    STRINTMAP* newEdgeStatus=new STRINTMAP;
    affectedOGPairs[ogpairKey]=newEdgeStatus;
    for(INTINTMAP_ITER cIter=currEdgeStatus->begin();cIter!=currEdgeStatus->end();cIter++)
    {
        (*newEdgeStatus)[speciesIDNameMap[cIter->first]]=cIter->second;
    }
    string& species=speciesIDNameMap[specID];   
    //(*newEdgeStatus)[species]=0;
    (*newEdgeStatus)[species]=(*newEdgeStatus)[species]-1;
    if((*newEdgeStatus)[species]>0)
    {
        cout <<"Found a duplicate for "<< species << " gene " << v->getName() << " ogkey " << ogpairKey << endl;
    }
    
    return 0;
}*/

// compute the regression weight during output
int
MetaLearner::dumpAllGraphs(int currK,int foldid)
{
    cout <<"MetaLearner::dumpAllGraphs" << endl;
    char foldoutDirName[1024];
    char aFName[1024];
    for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    {
        //cout <<"Dumping " << eIter->first << endl;
        SpeciesDataManager* sdm=eIter->second;
        const char* dirname=sdm->getOutputLoc();
        sprintf(foldoutDirName,"%s/fold%d",dirname,foldid);
        sprintf(aFName,"%s/var_mb_pw_k%d.txt",foldoutDirName,currK);
        ofstream oFile(aFName);
        FactorGraph* fg=sdm->getFactorGraph();
        VSET& varSet=sdm->getVariableManager()->getVariableSet();
        //fg->dumpVarMB_PairwiseFormat(oFile,varSet);
        map<int,SlimFactor*>& factorSet=fg->getAllFactors();
        PotentialManager* potMgr=sdm->getPotentialManager();
        for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
        {
            SlimFactor* sFactor=aIter->second;
            potMgr->dumpVarMB_PairwiseFormat(sFactor, oFile, varSet);
        }
        oFile.close();
        /*sprintf(aFName,"%s/net_ogspace_k%d.txt",foldoutDirName,currK);
        ofstream eFile(aFName);
        map<int,SlimFactor*>& factorSet=fg->getAllFactors();
        for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
        {
            SlimFactor* sFactor=aIter->second;
            INTDBLMAP& mbWts=sFactor->mbWts;
            for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
            {
                if(mbWts.find(mIter->first)==mbWts.end())
                {
                    cout <<"No neighbor " << mIter->first << " for " << sFactor->fId << endl;
                    exit(0);
                }
                double wtval=mbWts[mIter->first];
                int regulatorogid=ogr->getMappedOrthogroupID(varSet[mIter->first]->getName().c_str(),eIter->first.c_str());
                int targetogid=ogr->getMappedOrthogroupID(varSet[aIter->first]->getName().c_str(),eIter->first.c_str());
                //eFile << varSet[mIter->first]->getName()<< "\t" 
                //<< varSet[sFactor->vIds[0]]->getName() << "\t" << wtval 
                eFile << "\t" << regulatorogid << "\t" << targetogid << "\t" << wtval<< endl;
            }
        }
        eFile.close();*/
    }
    return 0;
}

int
MetaLearner::showModelParameters(int foldid)
{
    char foldoutDirName[1024];
    char aFName[1024];
    for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    {
        SpeciesDataManager* sdm=eIter->second;
        const char* dirname=sdm->getOutputLoc();
        sprintf(foldoutDirName,"%s/fold%d",dirname,foldid);
        sprintf(aFName,"%s/modelparams.txt",foldoutDirName);
        ofstream oFile(aFName);
        FactorGraph* fg=sdm->getFactorGraph();
        PotentialManager* potMgr=sdm->getPotentialManager();
        map<int,SlimFactor*>& slimFactorSet=fg->getAllFactors();
        VSET& varSet=sdm->getVariableManager()->getVariableSet();
        for(map<int,SlimFactor*>::iterator sIter=slimFactorSet.begin();sIter!=slimFactorSet.end();sIter++)
        {
            SlimFactor* sFactor=sIter->second;
            /*Potential* sPot=new Potential;
            sPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
            for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
            {
                Variable* aVar=varSet[mIter->first];
                sPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
            }
            sPot->potZeroInit();
            potMgr->populatePotential(sPot,false);
            sPot->initMBCovMean();
            double mbcondvar=sPot->getCondVariance();
            double mbbias=sPot->getCondBias();
            INTDBLMAP& mbwt=sPot->getCondWeight();*/
            //The weight of the model is really a space filler.
            double mbcondvar=0;
            double mbbias=0;
            unordered_map<int,double> mbwt;
            potMgr->computePotentialMBCovMean(sFactor, mbcondvar, mbbias, mbwt);
            Variable* var=varSet[sFactor->fId];
            oFile<<"Var="<<var->getName()<<"\tWt=-1"<<"\tCondVar="<<mbcondvar << "\tCondBias="<<mbbias<<"\tCondWt=";
            for(auto dIter=mbwt.begin();dIter!=mbwt.end();dIter++)
            {
                if(dIter!=mbwt.begin())
                {
                    oFile <<",";
                }
                Variable* mbVar=varSet[dIter->first];
                oFile<< mbVar->getName()<<"="<<dIter->second;
            }
            oFile << endl;
            mbwt.clear();
        }
        oFile.close();
    }
    return 0;
}



bool 
MetaLearner::checkMBSize(int cid, int u, int currK)
{
    bool check=true;
    SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[cid]];
    VSET& varSet=sdm->getVariableManager()->getVariableSet();
    Variable* sVar=varSet[u];
    FactorGraph* fg=sdm->getFactorGraph();
    SlimFactor* sFactor=fg->getFactorAt(u);
    map<string,int>& regulatorSet=sdm->getRegulators();
    if(sFactor->mergedMB.size()>=currK)
    {
        check=false;
    }
    return check;
}




int 
MetaLearner::generateData(int sampleCnt, int burnin,int foldid)
{
    map<int,ofstream*> newdataFiles;
    gsl_rng* rndgen=gsl_rng_alloc(gsl_rng_default);
    map<int,map<int,Potential*>*> potSet;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {
        int gid=speciesNameIDMap[gIter->first];
        FactorGraph* fg=gIter->second->getFactorGraph();
        VariableManager* varMgr=gIter->second->getVariableManager();
        VSET& varSet=varMgr->getVariableSet();
        map<int,Potential*>* pSet=new map<int,Potential*>;
        potSet[gid]=pSet;
        PotentialManager* potMgr=gIter->second->getPotentialManager();
        for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
        {
            int vId=vIter->second;
            SlimFactor* sFactor=fg->getFactorAt(vId);
            Potential* sPot=new Potential;
            (*pSet)[vId]=sPot;
            sPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
            for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
            {
                Variable* aVar=varSet[*mIter];
                sPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
            }
            sPot->potZeroInit();
            potMgr->populatePotential_Eff(sPot,false);
            sPot->initMBCovMean();
        }
    }
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        SpeciesDataManager* sdm=sIter->second;
        char fileName[1024];
        sprintf(fileName,"%s/newsamples_f%d.txt",sdm->getOutputLoc(),foldid);
        ofstream* oFile=new ofstream(fileName);
        newdataFiles[speciesNameIDMap[sIter->first]]=oFile;
        //Write the model file
        sprintf(fileName,"%s/newmodel_f%d.txt",sdm->getOutputLoc(),foldid);
        ofstream mFile(fileName);
        mFile <<"NodeCnt\t"<< subgraphVarSet.size() << endl;
        mFile <<"ContinuousNodes";
        int nodeId=0;
        for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
        {
            mFile <<"\t"<< nodeId;
            nodeId++;
        }
        mFile << endl;
        mFile <<"DiscreteNodes"<< endl;
        nodeId=0;
        for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
        {
            mFile <<"NodeName="<< vIter->first<< "\tNodeID="<< nodeId<< "\tParents=" 
            << "\tChildren=" <<"\tValues=0,1,2" << endl;
            nodeId++;
        }
        mFile.close();
    }

    for(map<int,map<int,Potential*>*>::iterator psIter=potSet.begin();psIter!=potSet.end();psIter++)
    {
        int did=psIter->first;
        INTDBLMAP initialSample;
        INTDBLMAP saveSample;
        int iter=0;
        map<int,Potential*>* pSet=psIter->second;
        generateInitSample(initialSample,rndgen,pSet);

        while(iter<(sampleCnt+burnin))
        {
            int selectFirst=0;
            //We just need to align localmodelid with the outputfile to which the data is generated
            if(iter>=burnin)
            {
                if(iter==burnin)
                {
                    cout <<"Starting to write"<< endl;
                }
                writeEvidenceTo(newdataFiles[did],initialSample);
            }
            for(map<string,int>::iterator vIter=subgraphVarSet.begin(); vIter!=subgraphVarSet.end();vIter++)
            {
                saveSample[vIter->second]=initialSample[vIter->second];
                FactorGraph* fg=speciesDataSet[speciesIDNameMap[psIter->first]]->getFactorGraph();
                SlimFactor* sFactor=fg->getFactorAt(vIter->second);
                Potential* sPot=(*pSet)[vIter->second];
                double sample_s=sPot->generateSample(initialSample,sFactor->fId,rndgen);
                int siter=0;
                while((sample_s<-160|| sample_s>160) && (siter<50))
                {
                    sample_s=sPot->generateSample(initialSample,sFactor->fId,rndgen);
                    siter++;
                }
                if(sample_s<-160)
                {
                    sample_s=-160;
                }
                else if(sample_s>160)
                {
                    sample_s=160;
                }
                if((isnan(sample_s))||(isinf(sample_s)))
                {
                    cout << "Found Nan/Inf for " << sFactor->fId <<"="<<sample_s  <<" Neighbors: ";
                    for(auto mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
                    {
                        cout << " "<< *mbIter <<"=" << initialSample[*mbIter];
                    }
                    cout << endl;
                }
                initialSample[sFactor->fId]=sample_s;
            }
            iter++;
        }
    }
    gsl_rng_free(rndgen);
    cout <<"Bottom clips " << endl;
    for(map<string,int>::iterator sIter=bottomClip.begin();sIter!=bottomClip.end();sIter++)
    {
        cout << sIter->first.c_str() << " " << sIter->second << endl;
    }
    cout <<"Top clips " << endl;
    for(map<string,int>::iterator sIter=topClip.begin();sIter!=topClip.end();sIter++)
    {
        cout << sIter->first.c_str() << " " << sIter->second << endl;
    }
    return 0;
}



int
MetaLearner::generateInitSample(INTDBLMAP& initialSample,gsl_rng* rndgen,map<int,Potential*>* pSet)
{
    double globalMean=0;
    globalVar=0;
    for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
    {
        int vId=vIter->second;
        Potential* sPot=(*pSet)[vId];
        cout <<vId <<"\t" << sPot->getCondBias() <<"\t" << sPot->getCondVariance() << endl;
        globalMean=globalMean+sPot->getCondBias();
        globalVar=globalVar+sPot->getCondVariance();
    }
    double norm=(double)(subgraphVarSet.size()*speciesDataSet.size());
    globalMean=globalMean/norm;
    globalVar=globalVar/norm;
    for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
    {
        double sampleval=gsl_ran_gaussian(rndgen,sqrt(globalVar));
        sampleval=sampleval+globalMean;
        initialSample[vIter->second]=sampleval;
    }
    return 0;
}


int
MetaLearner::writeEvidenceTo(ofstream* oFile,INTDBLMAP& sample)
{
    int nodeId=0;
    for(map<string,int>::iterator dIter=subgraphVarSet.begin();dIter!=subgraphVarSet.end();dIter++)
    {
        if(dIter!=subgraphVarSet.begin())
        {
            (*oFile)<<"\t";
        }
        double val=sample[dIter->second];
        if(val<-160)
        {
            val=-160;
            if(bottomClip.find(dIter->first)==bottomClip.end())
            {
                bottomClip[dIter->first]=1;
            }
            else
            {
                bottomClip[dIter->first]=bottomClip[dIter->first]+1;
            }
        }
        else if(val > 160)
        {
            val=160;
            if(topClip.find(dIter->first)==topClip.end())
            {
                topClip[dIter->first]=1;
            }
            else
            {
                topClip[dIter->first]=topClip[dIter->first]+1;
            }
        }
        (*oFile)<< nodeId << "=["<< exp(val) << "]";
        nodeId++;
    }
    (*oFile) << endl;
    return 0;
}
