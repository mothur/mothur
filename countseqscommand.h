#ifndef COuNTSEQSCOMMAND_H
#define COuNTSEQSCOMMAND_H

/*
 *  countseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 6/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "groupmap.h"

class CountSeqsCommand : public Command {
	
public:
	
	CountSeqsCommand(string);
	CountSeqsCommand();	
	~CountSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "count.seqs";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Count.seqs"; }
	string getDescription()		{ return "counts the number of sequences represented by each unique sequence in a namesfile"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
    
    struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};
    
	string namefile, groupfile, outputDir, groups;
	bool abort, large;
	vector<string> Groups, outputNames;
    int processors;
    
    int processSmall(string);
    int processLarge(string);
    map<int, string> processNameFile(string);
    map<int, string> getGroupNames(string, set<string>&);
    
    int createProcesses(GroupMap*&, string);
    int driver(unsigned long long, unsigned long long, string, GroupMap*&);
    
};

/***********************************************************************/
struct countData {
    unsigned long long start;
	unsigned long long end;
	MothurOut* m;
    string outputFileName, namefile, groupfile;
    GroupMap* groupMap;
    int total;
    vector<string> Groups;
    
	countData(){}
	countData(string fn, GroupMap* g, MothurOut* mout, unsigned long long st, unsigned long long en, string gfn, string nfn, vector<string> gr) {
        m = mout;
		start = st;
		end = en;
        groupMap = g;
        groupfile = gfn;
        namefile = nfn;
        outputFileName = fn;
        Groups = gr;
        total = 0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyCountThreadFunction(LPVOID lpParam){
	countData* pDataArray;
	pDataArray = (countData*)lpParam;
	try {
        ofstream out;
        pDataArray->m->openOutputFile(pDataArray->outputFileName, out);
        
        ifstream in;
		pDataArray->m->openInputFile(pDataArray->namefile, in);
		in.seekg(pDataArray->start);
        
        //print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			in.seekg(0);
		}else { //this accounts for the difference in line endings.
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in);
		}
        
        pDataArray->total = 0;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
            
			if (pDataArray->m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol; pDataArray->m->gobble(in); in >> secondCol; pDataArray->m->gobble(in);
            //cout << firstCol << '\t' << secondCol << endl;
            pDataArray->m->checkName(firstCol);
            pDataArray->m->checkName(secondCol);
            
			vector<string> names;
			pDataArray->m->splitAtChar(secondCol, names, ',');
			
			if (pDataArray->groupfile != "") {
				//set to 0
				map<string, int> groupCounts;
				int total = 0;
				for (int i = 0; i < pDataArray->Groups.size(); i++) { groupCounts[pDataArray->Groups[i]] = 0; }
				
				//get counts for each of the users groups
				for (int i = 0; i < names.size(); i++) {
					string group = pDataArray->groupMap->getGroup(names[i]);
					
					if (group == "not found") { pDataArray->m->mothurOut("[ERROR]: " + names[i] + " is not in your groupfile, please correct."); pDataArray->m->mothurOutEndLine(); }
					else {
						map<string, int>::iterator it = groupCounts.find(group);
						
						//if not found, then this sequence is not from a group we care about
						if (it != groupCounts.end()) {
							it->second++;
							total++;
						}
					}
				}
				
				if (total != 0) {
					out << firstCol << '\t' << total << '\t';
					for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
						out << it->second << '\t';
					}
					out << endl;
				}
			}else {
				out << firstCol << '\t' << names.size() << endl;
			}
			
			pDataArray->total += names.size();
		}
		in.close();
        out.close();

                
        return 0;
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "CountSeqsCommand", "MyCountThreadFunction");
		exit(1);
	}
}
#endif



#endif


