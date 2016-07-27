#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <map>
#include "CommonTypes.H"
using namespace std;


class VariableManager;
class EvidenceManager;
class PotentialManager;
class MetaMove;
class SpeciesDistance;
class MappedOrthogroupReader;
class MetaLearner
{
	public:
	MetaLearner();
	~MetaLearner();
	int setInputFName(const char*);
	int setMaxFactorSize(int);
	int setMaxFactorSize_Approx(int);
	int setMBCntForApproxGreedy(int);
	int setPenalty(double);
	int setConvergenceThreshold(double);
	int setRestrictedList(const char*);
	int setLambdaLearnSize(int);
	int setLambdaLearn();
	int setTrueGraph(const char* aGraphName);
	int setScoreOpt_Entropy();
	int setScoreOpt_PLL();
	int setInputOGList(const char*);

	int setOrthogroupReader(MappedOrthogroupReader*);
	int setSpeciesDistances(SpeciesDistance*);
	int setBeta1(double);
	int setBeta2(double);
	int init();
	int doCrossValidation(int);
	int start(int);
	int showModelParameters(int);
	int setPreRandomizeSplit();
	int setPreRandomizeSplitSeed(int seed);

	private:
		
	double getInitScore();
	double getScore();
	double getPriorDelta();
	int clearFoldSpecData();
	int initEdgeSet();
	double precomputeEmptyGraphPrior();
	int precomputePerSpeciesPrior();
	int initCondsetMap();
	int initCondsetMap_Nopool();
	INTINTMAP* getConditionSet(int cind);

	int collectMoves(int);
	int collectMoves_Orthogroups(int);
	int collectMoves_Deletions();
	int sortMoves();
	int makeMoves();

	int makeDelMoves();
	int attemptMove(MetaMove*,map<int,INTINTMAP*>&);
	int attemptDelMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars);
	int getNewPLLScore(int,INTINTMAP& conditionSet,Variable* u, Variable* v,double&,double&,int,INTDBLMAP&);
	double getPLLScore_Condition(string& condKey,SlimFactor* sFactor);
	double getPLLScore_Condition(string& condKey,SlimFactor* sFactor,int&,INTDBLMAP&);
	double getValidationPLLScore_Condition(int,int);
	double getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo);
	double getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo,int);
	double getEdgePrior_PerSpecies(int tfID, int targetID,SpeciesDataManager*);
	int genCondSetKey(INTINTMAP& condSet, string& aKey);
	int cleanUp();
	int dumpAllGraphs(int,int);
	bool checkMBSize(int,int,int);

	int generateData(int sampleCnt, int burnin,int foldid);
	int generateInitSample(INTDBLMAP& initialSample,gsl_rng* rndgen,map<int,Potential*>*);
	int writeEvidenceTo(ofstream* oFile,INTDBLMAP& sample);
	int initCondsetMap_Tree(SpeciesDistance::Species* node);
	int getLeaves(SpeciesDistance::Species* node, map<string,int>& leaves);

	map<int,string> speciesIDNameMap;
	map<string,int> speciesNameIDMap;
	map<string,SpeciesDataManager*> speciesDataSet;


	map<string,PotentialManager*> pooledPotentials;
	map<string,EvidenceManager*> pooledData;
	map<string,int> pooledKeysForDel;
	int maxFactorSize;
	int maxFactorSizeApprox;
	int maxMBCntGreedy;
	char inputFName[1024];
	double penalty;
	double convThreshold;
	char restrictedFName[1024];
	int lambdaLearnSize;
	bool learnLambda;

	map<string,INTINTMAP*> edgeConditionMap;
	map<string,STRINTMAP*> affectedOGPairs;
	map<string,double> ogpairPrior;
	FactorGraph* fg;
	vector<MetaMove*> moveSet;
	map<string,double> edgeValue;
	map<string,int> sharedEdges;
	char sharedEdgesFName[1024];
	map<int,INTINTMAP*> condsetMap;
	map<int,INTINTMAP*> condsetMap_Tree;
	map<string,int> condsetKeyIDMap;
	map<int,string> condsetIDKeyMap;
	map<int,int> restrictedOrthoList;
	bool dataPool;
	bool usePLL;
	bool useEntropy;
	map<int,INTDBLMAP*> validationPLLs;
	char trueGraphFName[1024];
	map<int,double> finalScores;
	map<string,int> subgraphVarSet;
	map<string,int> bottomClip;
	map<string,int> topClip;
	int totalFoldCnt;
	bool preRandomizeSplit;
	int preRandSeed;
	double globalVar;
	MappedOrthogroupReader* ogr;
	SpeciesDistance* speciesData;
	char orthoMapFName[1024];
	char speciesOrder[1024];
	map<int,int> inputOGList;
	map<int,int> inputRegulatorOGs;
	map<string,map<int,double>*> varNeighborhoodPrior_PerSpecies;//This stores the empty graph
	map<string,map<string,double>*> edgePresenceProb_PerSpecies;//This stores the edge probability
	double beta1;
	double beta2;
};
#endif
