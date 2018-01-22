/* 
 * File:   kruskalwalliscommand.h
 * Author: kiverson
 *
 * Created on June 26, 2012, 11:07 AM
 */

#ifndef KRUSKALWALLISCOMMAND_H
#define	KRUSKALWALLISCOMMAND_H

#include "command.hpp"
#include "inputdata.h"
#include "designmap.h"

class KruskalWallisCommand : public Command {
   
public:
	    
	KruskalWallisCommand(string);	
	KruskalWallisCommand();
	~KruskalWallisCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "kruskal.wallis";			}
	string getCommandCategory()		{ return "Hypothesis Testing";	}
    string getOutputPattern(string);	
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Kruskal.wallis"; }
	string getDescription()		{ return "Non-parametric method for testing whether samples originate from the same distribution."; }
    
    struct groupRank {
        string group;
        double value;
        double rank;        
    };
    
    int execute(); 
	void help() { m->mothurOut(getHelpString()); }
    void assignRank(vector<groupRank>&);
    void assignValue(vector<groupRank>&);
    
    
private:
    bool abort, allLines;
    string outputDir, sharedfile, designfile, mclass;
    vector<string> outputNames;
    set<string> labels;
    
    int process(vector<SharedRAbundVector*>&, DesignMap&, vector<string>);
};

#endif	/* KRUSKALWALLISCOMMAND_H */

