/* 
 * File:   kruskalwalliscommand.h
 * Author: kiverson
 *
 * Created on June 26, 2012, 11:07 AM
 */

#ifndef KRUSKALWALLISCOMMAND_H
#define	KRUSKALWALLISCOMMAND_H

#include "command.hpp"

class KruskalWallisCommand : public Command {

public:
	
	KruskalWallisCommand(string);	
	KruskalWallisCommand();
	~KruskalWallisCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "kruskalwallis";			}
	string getCommandCategory()		{ return "Hypothesis Testing";	}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/kruskalwallis"; }
	string getDescription()		{ return "Non-parametric method for testing whether samples originate from the same distribution."; }
    
    int execute(); 
	void help() { m->mothurOut(getHelpString()); }
    multimap<double,double> getRank(vector<KruskalWallisCommand::groupRank>);
    
private:
    string outputDir;
    vector<int> counts;
    vector<double> rankSums;
    vector<double> rankMeans;
    
};

#endif	/* KRUSKALWALLISCOMMAND_H */

