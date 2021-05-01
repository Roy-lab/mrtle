#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include <map>
#include <string>
#include <vector>

using namespace std;

class Error;
class Variable;
class VariableManager
{
	public:
		VariableManager();
		~VariableManager();
		Error::ErrorCode readVariables(const char*);
		int getVarID(const char*);
		bool isValid(int,int);
		map<int,Variable*>& getVariableSet();
		Variable* getVariableAt(int);

		Error::ErrorCode readVariablesFromTable(vector<string>& inputLines);
        Error::ErrorCode readVariablesFromData(string,int);
	private:
		map<string,int> varNameIDMap;
		map<int,Variable*> variableSet;
		int globalVarId;
};
#endif
