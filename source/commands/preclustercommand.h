#ifndef PRECLUSTERCOMMAND_H
#define PRECLUSTERCOMMAND_H


/*
 *  preclustercommand.h
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sequence.hpp"
#include "sequenceparser.h"
#include "sequencecountparser.h"
#include "alignment.hpp"
#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"
#include "filters.h"

/************************************************************/
struct seqPNode {
	int numIdentical;
	Sequence seq;
    string filteredSeq;
	string names;
	bool active;
	int diffs;
	seqPNode() {}
	seqPNode(int n, Sequence s, string nm) : numIdentical(n), seq(s), names(nm), active(1) { diffs = 0; filteredSeq = "";}
	~seqPNode() {}
};
/************************************************************/
inline bool comparePriorityTopDown(seqPNode first, seqPNode second) {  
    if (first.numIdentical > second.numIdentical) { return true;  }
    //else if (first.numIdentical == second.numIdentical) {
        //if (first.seq.getName() > second.seq.getName()) { return true; }
    //}

    return false; 
}
/************************************************************/
inline bool comparePriorityDownTop(seqPNode first, seqPNode second) {  
    if (first.numIdentical < second.numIdentical) { return true;  }
    //else if (first.numIdentical == second.numIdentical) {
        //if (first.seq.getName() > second.seq.getName()) { return true; }
    //}
    return false; 
}

//**********************************************************************************************************************
struct preClusterData {
    string fastafile;
    string namefile;
    string groupfile, countfile, method, align, newMName;
    OutputWriter* newFName;
    OutputWriter* newNName;
    MothurOut* m;
    int start;
    int end, count;
    int diffs, length;
    vector<string> groups;
    //vector<string> mapFileNames;
    bool topdown, hasCount, hasName;
    float match, misMatch, gapOpen, gapExtend;
    Utils util;
    vector<string> outputNames;
    map<string, vector<string> > outputTypes;
    vector<seqPNode> alignSeqs; //maps the number of identical seqs to a sequence
    Alignment* alignment;
    
    ~preClusterData() { if (alignment != NULL) { delete alignment; } }
    preClusterData(){}
    preClusterData(string f, string n, string g, string c, OutputWriter* nff,  OutputWriter* nnf, string nmf, vector<string> gr) {
        fastafile = f;
        groupfile = g;
        newFName = nff;
        newNName = nnf;
        newMName = nmf;
        groups = gr;
        hasName = false;
        namefile = n; if (namefile != "") { hasName = true; }
        hasCount = false;
        countfile = c; if (countfile != "") { hasCount = true; }
        count=0;
        m = MothurOut::getInstance();
    }
    void setVariables(int st, int en, int d, bool td, string me, string al, float ma, float misma, float gpOp, float gpEx) {
        start = st;
        end = en;
        diffs = d;
        topdown = td;
        method = me;
        align = al;
        match = ma;
        misMatch = misma;
        gapExtend = gpEx;
        gapOpen = gpOp;
        length = 0;
        
        if (method == "unaligned") {
            if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, 1000);	}
            else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 1000);			}
            else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
            else if(align == "noalign")		{	alignment = new NoAlign();													}
            else {
                m->mothurOut(align + " is not a valid alignment option. I will run the command using needleman.");
                m->mothurOutEndLine();
                alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 1000);
            }
        }else { alignment = NULL; }
    }
};

//************************************************************/
class PreClusterCommand : public Command {
	
public:
	PreClusterCommand(string);
	PreClusterCommand();
	~PreClusterCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pre.cluster";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Gevers D, Westcott SL (2011).  Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies.  PLoS ONE.  6:e27310.\nhttp://www.mothur.org/wiki/Pre.cluster"; }
	string getDescription()		{ return "implements a pseudo-single linkage algorithm with the goal of removing sequences that are likely due to pyrosequencing errors"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    CountTable ct;
	int diffs, length, processors;
    float match, misMatch, gapOpen, gapExtend;
	bool abort, bygroup, topdown;
	string fastafile, namefile, outputDir, groupfile, countfile, method, align;
	vector<string> outputNames;
	
	void createProcessesGroups(string, string, string);
    int mergeGroupCounts(string, string, string);
    void print(string, string, preClusterData*);
};


/**************************************************************************************************/

#endif


