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
    
    struct groupRank {
        sstring group;
        double value;
        double rank;        
    };
    
    int execute(); 
	void help() { m->mothurOut(getHelpString()); }
    void assignRank(vector<groupRank>);
    
private:
    string outputDir;
    vector<int> counts;
    vector<double> rankSums;
    vector<double> rankMeans;
    
  
        
    bool comparevalue(const groupRank &a, const groupRank &b) { return a.value < b.value; }
    bool equalvalue(const groupRank &a, const groupRank &b) { return a.value == b.value; }
    bool comparerank(const groupRank &a, const groupRank &b) { return a.rank < b.rank; }
    bool equalrank(const groupRank &a, const groupRank &b) { return a.rank == b.rank; }
    bool equalgroup(const groupRank &a, const groupRank &b) { return a.group == b.group; }
    
};

#endif	/* KRUSKALWALLISCOMMAND_H */

