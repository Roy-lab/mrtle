#include <fstream>
#include <iostream>
#include <cstring>

#include <unistd.h>

using namespace std;
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"

#include "Vertex.H"
#include "Graph.H"

#include "FactorGraph.H"
#include "PotentialManager.H"
#include "MetaMove.H"

#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "SpeciesDistance.H"
#include "SpeciesDataManager.H"
#include "MetaLearner.H"
#include "Framework.H"

Framework::Framework()
{
	epsThreshold=-1;
}

Framework::~Framework()
{
}

//We will use getopt here
//The options are 
//-m modelname
//-o outputdir
//-e epsilon to control the number of standard deviations above random
//-k maxfactorsize
//-s number of samples for approximate information estimation
//-x k for which we approximate information
//-n cnt of the top candidate MBs to save

Error::ErrorCode
Framework::init(int argc, char** argv)
{
	int optret='-';
	opterr=1;
	int oldoptind=optind;
	char orthoMapFName[1024];
	char speciesOrder[1024];
	while(optret=getopt(argc,argv,"f:k:x:p:t:v:l:o:g:r:d:m:s:n:b:q:")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return Error::UNKNOWN;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind ==2)
		{
			my_optarg=argv[oldoptind+1];	
		}
		else
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 'f':
			{
				metaLearner.setInputFName(my_optarg);
				break;
			}
			case 'k':
			{
				int aSize=atoi(my_optarg);
				metaLearner.setMaxFactorSize(aSize);
				break;
			}
			case 'x':
			{
				int aSize=atoi(my_optarg);
				metaLearner.setMaxFactorSize_Approx(aSize);
				break;
			}
			case 'p':
			{
				double penalty=atof(my_optarg);
				metaLearner.setPenalty(penalty);	
				break;
			}
			case 't':
			{
				double convThreshold=atof(my_optarg);
				metaLearner.setConvergenceThreshold(convThreshold);
				break;
			}
			case 'v':
			{
				cvCnt=atoi(my_optarg);
				break;
			}
			case 'l':
			{
				metaLearner.setRestrictedList(my_optarg);
				break;
			}
			case 'o':
			{
				if(strcmp(my_optarg,"entropy")==0)
				{
					metaLearner.setScoreOpt_Entropy();
				}
				else if(strcmp(my_optarg,"pll")==0)
				{
					metaLearner.setScoreOpt_PLL();
				}
				else
				{
					cout <<"Undefined score type "<< endl;
					return Error::UNKNOWN;
				}
				break;
			}
			case 'g':
			{
				metaLearner.setTrueGraph(my_optarg);
				break;
			}
			case 'r':
			{
				if(strcmp(my_optarg,"yes")==0)
				{
					metaLearner.setPreRandomizeSplit();
				}
				else if(isdigit(my_optarg[0]))
				{	
					metaLearner.setPreRandomizeSplit();
					metaLearner.setPreRandomizeSplitSeed(atoi(my_optarg));
				}
				break;
			}
			case 'd':
			{
				spData.readSpeciesTree(my_optarg);
				break;
			}
			case 'm':
			{
				strcpy(orthoMapFName,my_optarg);
				break;
			}
			case 's':
			{
				strcpy(speciesOrder,my_optarg);
				break;
			}
			case 'n':
			{
				metaLearner.setInputOGList(my_optarg);
				break;
			}
			case 'b':
			{
				double beta1 = atof(my_optarg);
				metaLearner.setBeta1(beta1);
				break;
			}
			case 'q':
			{
				double beta2 = atof(my_optarg);
				metaLearner.setBeta2(beta2);
				break;
			}
			default:
			{
				cout <<"Unhandled option " << c  << endl;
				return Error::UNKNOWN;
			}
		}
		oldoptind=optind;
	}
	metaLearner.init();
	mor.readSpeciesMapping(speciesOrder);
	mor.readFile(orthoMapFName);
	metaLearner.setOrthogroupReader(&mor);
	metaLearner.setSpeciesDistances(&spData);
	return Error::SUCCESS;
}

int 
Framework::start()
{
	//metaLearner.start();
	metaLearner.doCrossValidation(cvCnt);
        return 0;
}


int
main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout <<"factorGraphInf " <<  endl
			<<"-f modelname " << endl
			<< "-k maxfactorsize " << endl
			 << "-x maxfactorsize_approx" << endl
			 << "-p penalty" << endl
			 << "-t convergence_threshold" << endl
			 << "-v cross_validation_cnt" << endl
			 << "-l restrictedfname" << endl
			 << "-o optimizationcrit" <<  endl
			 << "-g knowngraph" << endl
			 << "-r randomseed" << endl
			 << "-d distancematrix" << endl
			 << "-m synergy_proper_orthogrps" << endl
			 << "-s tree" << endl
			<< "-n input_og_list" << endl
			<< "-b beta1" << endl
			<< "-q beta2" << endl;

		return 0;
	}
	Framework fw;
	if(fw.init(argc,argv)!=Error::SUCCESS)
	{
		return 0;
	}
	fw.start();
	return 0;
	
}

