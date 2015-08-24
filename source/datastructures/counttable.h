#ifndef Mothur_counttable_h
#define Mothur_counttable_h


//
//  counttable.h
//  Mothur
//
//  Created by Sarah Westcott on 6/26/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

//This class is designed to read a count table file and store its data.
//count table files look like:

/*
 Representative_Sequence	total	F003D000	F003D002	F003D004	F003D006	F003D008	F003D142	F003D144	F003D146	F003D148	F003D150	MOCK.GQY1XT001	
 GQY1XT001C296C	6051	409	985	923	937	342	707	458	439	387	464	0	
 GQY1XT001A3TJI	4801	396	170	413	442	306	769	581	576	497	651	0	
 GQY1XT001CS2B8	3018	263	226	328	460	361	336	248	290	187	319	0	
 GQY1XT001CD9IB	2736	239	177	256	405	306	286	263	248	164	392	0	
 
 or if no group info was used to create it
 
 Representative_Sequence	total	
 GQY1XT001C296C	6051
 GQY1XT001A3TJI	4801
 GQY1XT001CS2B8	3018
 GQY1XT001CD9IB	2736
 GQY1XT001ARCB1	2183
 GQY1XT001CNF2P	2796
 GQY1XT001CJMDA	1667
 GQY1XT001CBVJB	3758
 
 
 */


#include "mothurout.h"
#include "listvector.hpp"
#include "groupmap.h"

class CountTable {
    
    public:
    
        CountTable() { m = MothurOut::getInstance(); hasGroups = false; total = 0; uniques = 0; }
        ~CountTable() {}
    
        //reads and creates smart enough to eliminate groups with zero counts 
        int createTable(set<string>&, map<string, string>&, set<string>&); //seqNames, seqName->group, groupNames 
        int createTable(string, string, bool); //namefile, groupfile, createGroup
        int readTable(string, bool, bool);
    
        int printTable(string);
        int printHeaders(ofstream&);
        int printSeq(ofstream&, string);
        bool testGroups(string file); //used to check if file has group data without reading it.
        int copy(CountTable*);
    
        bool hasGroupInfo() { return hasGroups; }
        int getNumGroups() { return groups.size(); }
        vector<string> getNamesOfGroups() {  return groups;   }  //returns group names, if no group info vector is blank.
        int addGroup(string);
        int removeGroup(string);
        
        int renameSeq(string, string); //used to change name of sequence for use with trees
        int setAbund(string, string, int); //set abundance number of seqs for that group for that seq
        int push_back(string); //add a sequence 
        int push_back(string, int); //add a sequence 
        int push_back(string, vector<int>); //add a sequence with group info
        int remove(string); //remove seq
        int get(string); //returns unique sequence index for reading distance matrices like NameAssignment
        int size() { return indexNameMap.size(); }
    
        vector<string> getGroups(string); //returns vector of groups represented by this sequences
        vector<int> getGroupCounts(string);  //returns group counts for a seq passed in, if no group info is in file vector is blank. Order is the same as the groups returned by getGroups function.
        int getGroupCount(string, string); //returns number of seqs for that group for that seq
        int getGroupCount(string); // returns total seqs for that group
        int getNumSeqs(string); //returns total seqs for that seq, 0 if not found
        int setNumSeqs(string, int); //set total seqs for that seq, return -1 if not found
        int getNumSeqs() { return total; } //return total number of seqs
        int getNumUniqueSeqs() { return uniques; } //return number of unique/representative seqs
        int getGroupIndex(string); //returns index in getGroupCounts vector of specific group
    
        vector<string> getNamesOfSeqs();
        vector<string> getNamesOfSeqs(string);
        int mergeCounts(string, string); //combines counts for 2 seqs, saving under the first name passed in.
        ListVector getListVector();
        map<string, int> getNameMap();
    
    private:
        string filename;
        MothurOut* m;
        bool hasGroups;
        int total, uniques;
        vector<string> groups;
        vector< vector<int> > counts;
        vector<int> totals;
        vector<int> totalGroups;
        map<string, int> indexNameMap;
        map<string, int> indexGroupMap;
    
};

#endif
