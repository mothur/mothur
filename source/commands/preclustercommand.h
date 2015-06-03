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


/************************************************************/
struct seqPNode {
	int numIdentical;
	Sequence seq;
	string names;
	bool active;
	int diffs;
	seqPNode() {}
	seqPNode(int n, Sequence s, string nm) : numIdentical(n), seq(s), names(nm), active(1) { diffs = 0; }
	~seqPNode() {}
};
/************************************************************/
inline bool comparePriorityTopDown(seqPNode first, seqPNode second) {  
    if (first.numIdentical > second.numIdentical) { return true;  }
    else if (first.numIdentical == second.numIdentical) { 
        if (first.seq.getName() > second.seq.getName()) { return true; }
    }
    return false; 
}
/************************************************************/
inline bool comparePriorityDownTop(seqPNode first, seqPNode second) {  
    if (first.numIdentical < second.numIdentical) { return true;  }
    else if (first.numIdentical == second.numIdentical) { 
        if (first.seq.getName() > second.seq.getName()) { return true; }
    }
    return false; 
}
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
    SequenceParser* parser;
    SequenceCountParser* cparser;
    CountTable ct;
    Alignment* alignment;
    
	int diffs, length, processors;
    float match, misMatch, gapOpen, gapExtend;
	bool abort, bygroup, topdown;
	string fastafile, namefile, outputDir, groupfile, countfile, method, align;
	vector<seqPNode> alignSeqs; //maps the number of identical seqs to a sequence
	map<string, string> names; //represents the names file first column maps to second column
	map<string, int> sizes;  //this map a seq name to the number of identical seqs in the names file
	map<string, int>::iterator itSize; 
//	map<string, bool> active; //maps sequence name to whether it has already been merged or not.
	vector<string> outputNames;
	
	int readFASTA();
	void readNameFile();
	//int readNamesFASTA();
	int calcMisMatches(string, string);
	void printData(string, string, string); //fasta filename, names file name
	int process(string);
	int loadSeqs(map<string, string>&, vector<Sequence>&, string);
	int driverGroups(string, string, string, int, int, vector<string> groups);
	int createProcessesGroups(string, string, string, vector<string>);
    int mergeGroupCounts(string, string, string);
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct preClusterData {
	string fastafile; 
	string namefile; 
	string groupfile, countfile;
	string newFName, newNName, newMName, method, align;
	MothurOut* m;
	int start;
	int end, count;
	int diffs, threadID;
	vector<string> groups;
	vector<string> mapFileNames;
    bool topdown;
	float match, misMatch, gapOpen, gapExtend;
    
	preClusterData(){}
	preClusterData(string f, string n, string g, string c, string nff,  string nnf, string nmf, vector<string> gr, MothurOut* mout, int st, int en, int d, bool td, int tid, string me, string al, float ma, float misma, float gpOp, float gpEx) {
		fastafile = f;
		namefile = n;
		groupfile = g;
		newFName = nff;
		newNName = nnf;
		newMName = nmf;
		m = mout;
		start = st;
		end = en;
		diffs = d;
		threadID = tid;
		groups = gr;
        countfile = c;
        topdown = td;
        count=0;
        method = me;
        align = al;
        match = ma;
        misMatch = misma;
        gapExtend = gpEx;
        gapOpen = gpOp;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyPreclusterThreadFunction(LPVOID lpParam){ 
	preClusterData* pDataArray;
	pDataArray = (preClusterData*)lpParam;
	
	try {
        
        Alignment* alignment;
        
        if(pDataArray->align == "gotoh")			{	alignment = new GotohOverlap(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch, 1000);	}
        else if(pDataArray->align == "needleman")	{	alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, 1000);			}
        else if(pDataArray->align == "blast")		{	alignment = new BlastAlignment(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch);		}
        else if(pDataArray->align == "noalign")		{	alignment = new NoAlign();													}
        else {
            pDataArray->m->mothurOut(align + " is not a valid alignment option. I will run the command using needleman.");
            pDataArray->m->mothurOutEndLine();
            alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, 1000);
        }
        
		//parse fasta and name file by group
		SequenceParser* parser;
        SequenceCountParser* cparser;
        if (pDataArray->countfile != "") {
            cparser = new SequenceCountParser(pDataArray->countfile, pDataArray->fastafile);
        }else {
            if (pDataArray->namefile != "") { parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile, pDataArray->namefile);	}
            else				{ parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile);			}
        }
        
 		int numSeqs = 0;
		vector<seqPNode> alignSeqs;
		//clear out old files
		ofstream outF; pDataArray->m->openOutputFile(pDataArray->newFName, outF); outF.close();
		ofstream outN; pDataArray->m->openOutputFile(pDataArray->newNName, outN);  outN.close();
		
		//precluster each group
		for (int k = pDataArray->start; k < pDataArray->end; k++) {
			
            pDataArray->count++;
            
			int start = time(NULL);
			
            if (pDataArray->m->control_pressed) {  delete parser; delete alignment;return 0; }
			
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("Processing group " + pDataArray->groups[k] + ":"); pDataArray->m->mothurOutEndLine();
			
			map<string, string> thisNameMap;
            vector<Sequence> thisSeqs;
			if (pDataArray->groupfile != "") { 
                thisSeqs = parser->getSeqs(pDataArray->groups[k]);
            }else if (pDataArray->countfile != "") {
                thisSeqs = cparser->getSeqs(pDataArray->groups[k]);
            }
			if (pDataArray->namefile != "") {  thisNameMap = parser->getNameMap(pDataArray->groups[k]); }
			
			//fill alignSeqs with this groups info.
			////////////////////////////////////////////////////
			//numSeqs = loadSeqs(thisNameMap, thisSeqs); same function below
			
			int length = 0;
            set<int> lengths;
			alignSeqs.clear();
			map<string, string>::iterator it;
			bool error = false;
            map<string, int> thisCount;
            if (pDataArray->countfile != "") { thisCount = cparser->getCountTable(pDataArray->groups[k]);  }

		 	
			for (int i = 0; i < thisSeqs.size(); i++) {
				
				if (pDataArray->m->control_pressed) { delete parser; delete alignment; return 0; }
				
				if (pDataArray->namefile != "") {
					it = thisNameMap.find(thisSeqs[i].getName());
					
					//should never be true since parser checks for this
					if (it == thisNameMap.end()) { pDataArray->m->mothurOut(thisSeqs[i].getName() + " is not in your names file, please correct."); pDataArray->m->mothurOutEndLine(); error = true; }
					else{
						//get number of reps
						int numReps = 1;
						for(int j=0;j<(it->second).length();j++){
							if((it->second)[j] == ','){	numReps++;	}
						}
						
						seqPNode tempNode(numReps, thisSeqs[i], it->second);
						alignSeqs.push_back(tempNode);
						lengths.insert(thisSeqs[i].getAligned().length());
					}	
				}else { //no names file, you are identical to yourself 
					int numRep = 1;
                    if (pDataArray->countfile != "") { 
                        map<string, int>::iterator it2 = thisCount.find(thisSeqs[i].getName());
                        
                        //should never be true since parser checks for this
                        if (it2 == thisCount.end()) { pDataArray->m->mothurOut(thisSeqs[i].getName() + " is not in your count file, please correct."); pDataArray->m->mothurOutEndLine(); error = true; }
                        else { numRep = it2->second;  }
                    }
                    seqPNode tempNode(numRep, thisSeqs[i], thisSeqs[i].getName());
                    alignSeqs.push_back(tempNode);
					lengths.insert(thisSeqs[i].getAligned().length());
				}
			}
			
            if (lengths.size() > 1) { pDataArray->method = "unaligned"; }
            else if (lengths.size() == 1) {  pDataArray->method = "aligned"; }
            
            length = *(lengths.begin());
			//sanity check
			if (error) { pDataArray->m->control_pressed = true; }
			
			thisSeqs.clear();
			numSeqs = alignSeqs.size();
			
			////////////////////////////////////////////////////
			
			if (pDataArray->m->control_pressed) {   delete parser; delete alignment; return 0; }
			
            if (pDataArray->method == "aligned") { if (pDataArray->diffs > length) { pDataArray->m->mothurOut("Error: diffs is greater than your sequence length."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true; delete alignment; return 0;  } }
			
			////////////////////////////////////////////////////
			//int count = process(); - same function below
			
			ofstream out;
			pDataArray->m->openOutputFile(pDataArray->newMName+pDataArray->groups[k]+".map", out);
			pDataArray->mapFileNames.push_back(pDataArray->newMName+pDataArray->groups[k]+".map");
			
            //sort seqs by number of identical seqs
            if (pDataArray->topdown) { sort(alignSeqs.begin(), alignSeqs.end(), comparePriorityTopDown);  }
            else {  sort(alignSeqs.begin(), alignSeqs.end(), comparePriorityDownTop);  }
            
			int count = 0;
			
            if (pDataArray->topdown) {
                //think about running through twice...
                for (int i = 0; i < numSeqs; i++) {
                    
                    //are you active
                    //			itActive = active.find(alignSeqs[i].seq.getName());
                    
                    if (alignSeqs[i].active) {  //this sequence has not been merged yet
                        
                        string chunk = alignSeqs[i].seq.getName() + "\t" + toString(alignSeqs[i].numIdentical) + "\t" + toString(0) + "\t" + alignSeqs[i].seq.getAligned() + "\n";

                        //try to merge it with all smaller seqs
                        for (int j = i+1; j < numSeqs; j++) {
                            
                            if (pDataArray->m->control_pressed) { delete parser; delete alignment; return 0; }
                            
                            if (alignSeqs[j].active) {  //this sequence has not been merged yet
                                //are you within "diff" bases
                                //int mismatch = calcMisMatches(alignSeqs[i].seq.getAligned(), alignSeqs[j].seq.getAligned());
                                ////////////////////////////////////////////////////
                                int mismatch = 0;
                                
                                if (pDataArray->method == "unaligned") {
                                    //align to eachother
                                    Sequence seqI("seq1", alignSeqs[i].seq.getAligned());
                                    Sequence seqJ("seq2", alignSeqs[j].seq.getAligned());
                                    
                                    //align seq2 to seq1 - less abundant to more abundant
                                    alignment->align(seqJ.getUnaligned(), seqI.getUnaligned());
                                    string seq2 = alignment->getSeqAAln();
                                    string seq1 = alignment->getSeqBAln();
                                    
                                    //chop gap ends
                                    int startPos = 0;
                                    int endPos = seq2.length()-1;
                                    for (int i = 0; i < seq2.length(); i++) {  if (isalpha(seq2[i])) { startPos = i; break; } }
                                    for (int i = seq2.length()-1; i >= 0; i--) {  if (isalpha(seq2[i])) { endPos = i; break; } }
                                    
                                    //count number of diffs
                                    for (int i = startPos; i <= endPos; i++) {
                                        if (seq2[i] != seq1[i]) { mismatch++; }
                                        if (mismatch > pDataArray->diffs) { mismatch = length; break;  } //to far to cluster
                                    }
                                }else {
                                    for (int k = 0; k < alignSeqs[i].seq.getAligned().length(); k++) {
                                        //do they match
                                        if (alignSeqs[i].seq.getAligned()[k] != alignSeqs[j].seq.getAligned()[k]) { mismatch++; }
                                        if (mismatch > pDataArray->diffs) { mismatch = length; break; } //to far to cluster
                                    }
                                }
                                ////////////////////////////////////////////////////
                                
                                if (mismatch <= pDataArray->diffs) {
                                    //merge
                                    alignSeqs[i].names += ',' + alignSeqs[j].names;
                                    alignSeqs[i].numIdentical += alignSeqs[j].numIdentical;
                                    
                                    alignSeqs[j].active = 0;
                                    alignSeqs[j].numIdentical = 0;
                                    alignSeqs[j].diffs = mismatch;
                                    count++;
                                    chunk += alignSeqs[j].seq.getName() + "\t" + toString(alignSeqs[j].numIdentical) + "\t" + toString(mismatch) + "\t" + alignSeqs[j].seq.getAligned() + "\n";
                                }
                            }//end if j active
                        }//end for loop j
                        
                        //remove from active list 
                        alignSeqs[i].active = 0;
                        
                        out << "ideal_seq_" << (i+1) << '\t' << alignSeqs[i].numIdentical << endl << chunk << endl;
                        
                    }//end if active i
                    if(i % 100 == 0)	{ pDataArray->m->mothurOutJustToScreen(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
                }
                
            }else {
                map<int, string> mapFile;
                map<int, int> originalCount;
                map<int, int>::iterator itCount;
                for (int i = 0; i < numSeqs; i++) { mapFile[i] = ""; originalCount[i] = alignSeqs[i].numIdentical; }
                
                //think about running through twice...
                for (int i = 0; i < numSeqs; i++) {
                    
                    //try to merge it into larger seqs
                    for (int j = i+1; j < numSeqs; j++) {
                        
                        if (pDataArray->m->control_pressed) { out.close(); delete alignment; return 0; }
                        
                        if (originalCount[j] > originalCount[i]) {  //this sequence is more abundant than I am
                            //are you within "diff" bases
                            //int mismatch = calcMisMatches(alignSeqs[i].seq.getAligned(), alignSeqs[j].seq.getAligned());
                            int mismatch = 0;
                            
                            if (pDataArray->method == "unaligned") {
                                //align to eachother
                                Sequence seqI("seq1", alignSeqs[i].seq.getAligned());
                                Sequence seqJ("seq2", alignSeqs[j].seq.getAligned());
                                
                                //align seq2 to seq1 - less abundant to more abundant
                                alignment->align(seqI.getUnaligned(), seqJ.getUnaligned());
                                string seq2 = alignment->getSeqAAln();
                                string seq1 = alignment->getSeqBAln();
                                
                                //chop gap ends
                                int startPos = 0;
                                int endPos = seq2.length()-1;
                                for (int i = 0; i < seq2.length(); i++) {  if (isalpha(seq2[i])) { startPos = i; break; } }
                                for (int i = seq2.length()-1; i >= 0; i--) {  if (isalpha(seq2[i])) { endPos = i; break; } }
                                
                                //count number of diffs
                                for (int i = startPos; i <= endPos; i++) {
                                    if (seq2[i] != seq1[i]) { mismatch++; }
                                    if (mismatch > pDataArray->diffs) { mismatch = length; break;  } //to far to cluster
                                }
                            }else {

                                for (int k = 0; k < alignSeqs[i].seq.getAligned().length(); k++) {
                                    //do they match
                                    if (alignSeqs[i].seq.getAligned()[k] != alignSeqs[j].seq.getAligned()[k]) { mismatch++; }
                                    if (mismatch > pDataArray->diffs) { mismatch = length; break; } //to far to cluster
                                }
                            }
                            if (mismatch <= pDataArray->diffs) {
                                //merge
                                alignSeqs[j].names += ',' + alignSeqs[i].names;
                                alignSeqs[j].numIdentical += alignSeqs[i].numIdentical;
                                
                                mapFile[j] = alignSeqs[i].seq.getName() + "\t" + toString(alignSeqs[i].numIdentical) + "\t" + toString(mismatch) + "\t" + alignSeqs[i].seq.getAligned() + "\n" + mapFile[i];
                                alignSeqs[i].numIdentical = 0;
                                originalCount.erase(i);
                                mapFile[i] = "";
                                count++;
                                j+=numSeqs; //exit search, we merged this one in.
                            }
                        }//end abundance check
                    }//end for loop j
                    
                    if(i % 100 == 0)	{ pDataArray->m->mothurOutJustToScreen(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
                }
                
                for (int i = 0; i < numSeqs; i++) {
                    if (alignSeqs[i].numIdentical != 0) {
                        out << "ideal_seq_" << (i+1) << '\t' << alignSeqs[i].numIdentical << endl  << alignSeqs[i].seq.getName() + "\t" + toString(alignSeqs[i].numIdentical) + "\t" + toString(0) + "\t" + alignSeqs[i].seq.getAligned() + "\n" << mapFile[i] << endl;
                    }
                }

            }
            out.close();
            if(numSeqs % 100 != 0)	{ pDataArray->m->mothurOut(toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count)); pDataArray->m->mothurOutEndLine();	}
			////////////////////////////////////////////////////
			
			if (pDataArray->m->control_pressed) {  delete parser; return 0; }
			
			pDataArray->m->mothurOut("Total number of sequences before pre.cluster was " + toString(alignSeqs.size()) + ".");pDataArray-> m->mothurOutEndLine();
			pDataArray->m->mothurOut("pre.cluster removed " + toString(count) + " sequences."); pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOutEndLine(); 
			
			////////////////////////////////////////////////////
			//printData(pDataArray->newFFile, pDataArray->newNFile); - same as below
			ofstream outFasta;
			ofstream outNames;
			
			pDataArray->m->openOutputFileAppend(pDataArray->newFName, outFasta);
			pDataArray->m->openOutputFileAppend(pDataArray->newNName, outNames);
						
			for (int i = 0; i < alignSeqs.size(); i++) {
				if (alignSeqs[i].numIdentical != 0) {
					alignSeqs[i].seq.printSequence(outFasta); 
					if (pDataArray->countfile != "") {  outNames << pDataArray->groups[k] << '\t' << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl; 
                    }else {  outNames << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl;  }

				}
			}
			
			outFasta.close();
			outNames.close();
			////////////////////////////////////////////////////
			
			pDataArray->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); pDataArray->m->mothurOutEndLine(); 
			
		}
		
        delete alignment;
        
		return numSeqs;
		

	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "PreClusterCommand", "MyPreclusterThreadFunction");
		exit(1);
	}
} 
#endif

/**************************************************************************************************/


#endif


