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
#include "sequence.hpp"
#include "sharedrabundvectors.hpp"

class CountTable {

    public:

        CountTable() { m = MothurOut::getInstance(); hasGroups = false; total = 0; uniques = 0; isCompressed = false; }
        ~CountTable() {}

        //reads and creates smart enough to eliminate groups with zero counts
        int createTable(set<string>&, map<string, string>&, set<string>&); //seqNames, seqName->group, groupNames
        int createTable(string, string, vector<string>, bool createGroup=false); //namefile, groupfile, selectedGroups, createGroup,
        int readTable(string, bool, bool); //filename, readGroups, mothurRunning
        int readTable(string, bool, bool, vector<string>); //filename, readGroups, mothurRunning, groups to save (if none provided, read all groups)
        int readTable(string, bool, bool, set<string>); //filename, readGroups, mothurRunning, namesofSeqs to save (if none provided, read all seqs)
        int readTable(string, string); //filename, format - if format=fasta, read fasta file and create unique table
    
        int zeroOutTable(); //turn all counts to zeros
        void eliminateZeroSeqs();
        int clearTable();
        bool isCountTable(string);
        bool isTableCompressed() { return isCompressed; }
        int copy(CountTable*); //copy countable
        bool inTable(string);  //accepts sequence name and returns true if sequence is in table, false if not present


        //all print commands ignore zeroed out seqs
        vector<string> printCompressedTable(string, vector<string> optionalGroups=nullVector); //nameOfFile, optionalVectorOfGroups (if empty, prints all possible groups), returns names of seqs in table - excludes zeroed reads
        vector<string> printTable(string); //preserves order in original, defaults compress to state of original file
        vector<string> printTable(string, bool compress); //preserves order in original, printing compressed or not based on compress flag pasted in
        vector<string> printSortedTable(string); //sorted by seqName
        int printHeaders(ofstream&, vector<string> optionalGroups=nullVector);
        vector<string> getHardCodedHeaders(); //Representative_Sequence, total
        int printSeq(ofstream&, string);
    
        bool testGroups(string file); //used to check if file has group data without reading it
        bool testGroups(string file, vector<string>&); //used to check if file has group data without reading it, return groups if found.
        bool hasGroupInfo() { return hasGroups; }
        int getNumGroups() { return (int)groups.size(); }
        vector<string> getNamesOfGroups() {  return groups;   }  //returns group names, if no group info vector is blank.
        bool setNamesOfGroups(vector<string>);
        int addGroup(string);
        int removeGroup(string); //pass in group name
        int removeGroup(int minSize);  //removes any groups with numSeqs < minSize

        int renameSeq(string, string); //used to change name of sequence for use with trees
        int setAbund(string, string, int); //set abundance number of seqs for that group for that seq
        int mergeCounts(string, string); //combines counts for 2 seqs, saving under the first name passed in.
        int push_back(string); //add a sequence
        int push_back(string, int); //add a sequence
        int push_back(string, vector<int>); //add a sequence with group info
        int push_back(string, vector<int>, bool); //add a sequence with group info, no error - ignore dups
        int remove(string); //remove seq
        int get(string); //returns unique sequence index for reading distance matrices like NameAssignment
        int size() { return (int)indexNameMap.size(); }

        vector<string> getGroups(string); //returns vector of groups represented by this sequence
        vector<int> getGroupCounts(string);  //returns group counts for a seq passed in, if no group info is in file vector is blank. Order is the same as the groups returned by getGroups function.
        int getGroupCount(string, string); //returns number of seqs for that group for that seq
        int getGroupCount(string); // returns total seqs for that group
        int getNumSeqs(string); //returns total seqs for that seq, 0 if not found
        int setNumSeqs(string, int); //set total seqs for that seq, return -1 if not found
        int getNumSeqs() { return total; } //return total number of seqs
        int getNumUniqueSeqs() { return uniques; } //return number of unique/representative seqs
        int getNumSeqsSmallestGroup(); //returns size of smallest group. If no groups, returns total num seqs (includes non uniques)

        vector<string> getNamesOfSeqs(); //return names of all seqeunce in table
        vector<string> getNamesOfSeqs(string); //returns names of seqs in specific group in table
        vector<string> getNamesOfSeqs(vector<string>); //returns names of seqs in specific set of groups in table
    
        ListVector getListVector();
        SharedRAbundVectors* getShared(map<string, string>&);
        SharedRAbundVectors* getShared(vector<string>, map<string, string>&); //set of groups selected
        map<string, int> getNameMap();  //sequenceName -> total number of sequences it represents
        map<string, int> getNameMap(string);  //sequenceName -> total number of sequences it represents in that group
    

    private:
        string filename;
        MothurOut* m;
        Utils util;
        bool hasGroups, isCompressed;
        int total, uniques;
        vector<string> groups;
        vector< vector<countTableItem> > counts; //countTableItem ((int)abund, (int)group). each line in counts represents a sequence line from the count table file(sparse). The vector<ountTableItem> are sorted by group, so that you can stop search early if group is not found. For example:  seq1 10 5 0 0 1 0 0 0 3 0 0 1 0 0 - 13 groups, but seq1 is only present in 4 samples. Let's save space by not storing 0 abunds. seq1's vector<ountTableItem> (5,0),(1,3),(3,7),(1,10). Group0 = 5, Group3 = 1, Group7 = 3, Group10 = 1.  
        vector<int> totals;
        vector<int> totalGroups;
        map<string, int> indexNameMap; //maps seqName -> vector index in counts. seq1 -> 1 would mean seq1's counts are stored in counts[1].
        map<string, int> indexGroupMap;
    
        int find(int seq, int group); //returns index of countTableItem for group passed in. If group is not present in seq, returns -1
        int getAbund(int seq, int group); //returns abundance of countTableItem for seq and group passed in. If group is not present in seq, returns 0
        vector<countTableItem> getItems(string); //returns group counts for a seq passed in, if no group info is in file vector is blank. sorted by group
        vector<int> expandAbunds(int index);
        vector<int> expandAbunds(vector<countTableItem>& items);
        vector<countTableItem> compressAbunds(vector<int> abunds);
        int printGroupAbunds(ofstream& out, int index);
        int sortCountTable();
        int sortRow(int);

};

#endif
