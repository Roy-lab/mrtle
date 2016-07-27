#ifndef _SPECIES_DATA_
#define _SPECIES_DATA_

#include <map>
#include <string>
using namespace std;
class SpeciesDistance
{
	public:
	SpeciesDistance();
	~SpeciesDistance();
	int readSpeciesTree(const char*);
	double getEdgeStatusProb(map<string,int>& edgeStatus);
	double getProbGainGain(string&);
	double getProbGainLoss(string&);
	double getProbLossGain(string&);
	double getProbLossLoss(string&);

	struct Species;

	struct Species 
	{
		string name;
		Species* parent;
		Species* leftchild;
		Species* rightchild;
		double p_maintain_edge;
		double p_gain;	
	};
	
	Species* getRoot();
	double getSubTree(bool,Species*,map<string,int>&);
	private:
	map<string,Species*> speciesSet;
	map<int,double> binEdgeKeyProbMap;
	Species* root;
};
#endif
