#ifndef CHOPSEQSCOMMAND_H
#define CHOPSEQSCOMMAND_H

/*
 *  chopseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 5/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sequence.hpp"
#include "qualityscores.h"

class ChopSeqsCommand : public Command {
	
	public:
	
		ChopSeqsCommand(string);
		ChopSeqsCommand();	
		~ChopSeqsCommand(){};
	
		vector<string> setParameters();
		string getCommandName()			{ return "chop.seqs";		}
		string getCommandCategory()		{ return "Sequence Processing"; }
		
        string getHelpString();	
        string getOutputPattern(string);	
		string getCitation() { return "http://www.mothur.org/wiki/Chops.seqs"; }
		string getDescription()		{ return "trim sequence length"; }
	
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }		
	
	private:
    
		string fastafile, outputDir, keep, namefile, groupfile, countfile, qualfile;
		bool abort, countGaps, Short, keepN;
		int numbases, processors;
		vector<string> outputNames;
		
		string getChopped(Sequence, string&);
        bool driver (linePair, string, string, string, string);
        bool createProcesses(vector<linePair>, string, string, string, string);
        int processQual(string, string);
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct chopData {
	string filename; 
	string outFasta, outAccnos, keep, qualfile, fastaFileTemp;
	unsigned long long start;
	unsigned long long end;
	int numbases, count;
    bool countGaps, Short, wroteAccnos, keepN;
	MothurOut* m;
	string namefile;
	map<string, int> nameMap;
	
	
	chopData(){}
	chopData(string f, string ff, string a, MothurOut* mout, unsigned long long st, unsigned long long en, string k, bool cGaps, int nbases, bool S, bool kn, string qu, string ft) {
		filename = f;
		outFasta = ff;
        outAccnos = a;
		m = mout;
		start = st;
		end = en;
        keep = k;
        countGaps = cGaps;
        numbases = nbases;
        Short = S;
		wroteAccnos = false;
        keepN = kn;
        qualfile = qu;
        fastaFileTemp = ft;
        
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyChopThreadFunction(LPVOID lpParam){ 
	chopData* pDataArray;
	pDataArray = (chopData*)lpParam;
	
	try {
        ofstream out;
		pDataArray->m->openOutputFile(pDataArray->outFasta, out);
        
        ofstream outAcc;
		pDataArray->m->openOutputFile(pDataArray->outAccnos, outAcc);
        
		ifstream in;
		pDataArray->m->openInputFile(pDataArray->filename, in);
        
        ofstream outfTemp;
        if (pDataArray->fastaFileTemp != "") { pDataArray->m->openOutputFile(pDataArray->fastaFileTemp, outfTemp); }

        
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			in.seekg(0);
            pDataArray->m->zapGremlins(in);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}

		bool done = false;
        bool wroteAccnos = false;
		pDataArray->count = 0;

		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
						
			if (pDataArray->m->control_pressed) {  in.close(); out.close(); outAcc.close(); pDataArray->m->mothurRemove(pDataArray->outFasta); pDataArray->m->mothurRemove(pDataArray->outAccnos); if (pDataArray->fastaFileTemp != "") { outfTemp.close(); pDataArray->m->mothurRemove(pDataArray->fastaFileTemp); } return 0;  }
            
            Sequence seq(in); pDataArray->m->gobble(in);
			
			if (seq.getName() != "") {
				//string newSeqString = getChopped(seq);
                ///////////////////////////////////////////////////////////////////////
                string qualValues = "";
                string temp = seq.getAligned();
                string tempUnaligned = seq.getUnaligned();
                
                if (pDataArray->countGaps) {
                    //if needed trim sequence
                    if (pDataArray->keep == "front") {//you want to keep the beginning
                        int tempLength = temp.length();
                        
                        if (tempLength > pDataArray->numbases) { //you have enough bases to remove some
                            
                            int stopSpot = 0;
                            int numBasesCounted = 0;
                            
                            for (int i = 0; i < temp.length(); i++) {
                                //eliminate N's
                                if (!pDataArray->keepN) { if (toupper(temp[i]) == 'N') { temp[i] = '.'; } }
                                
                                numBasesCounted++; 
                                
                                if (numBasesCounted >= pDataArray->numbases) { stopSpot = i; break; }
                            }
                            
                            if (stopSpot == 0) { temp = ""; }
                            else {  temp = temp.substr(0, stopSpot+1);  }
							
                        }else { 
                            if (!pDataArray->Short) { temp = ""; } //sequence too short
                        }
                    }else { //you are keeping the back
                        int tempLength = temp.length();
                        if (tempLength > pDataArray->numbases) { //you have enough bases to remove some
                            
                            int stopSpot = 0;
                            int numBasesCounted = 0;
                            
                            for (int i = (temp.length()-1); i >= 0; i--) {
                                //eliminate N's
                                if (!pDataArray->keepN) { if (toupper(temp[i]) == 'N') { temp[i] = '.'; } }
                                
                                numBasesCounted++; 
                                
                                if (numBasesCounted >= pDataArray->numbases) { stopSpot = i; break; }
                            }
                            
                            if (stopSpot == 0) { temp = ""; }
                            else {  temp = temp.substr(stopSpot+1);  }
                        }else { 
                            if (!pDataArray->Short) { temp = ""; } //sequence too short
                        }
                    }
                    
                }else{
                    
                    //if needed trim sequence
                    if (pDataArray->keep == "front") {//you want to keep the beginning
                        int tempLength = tempUnaligned.length();
                        
                        if (tempLength > pDataArray->numbases) { //you have enough bases to remove some
                            
                            int stopSpot = 0;
                            int numBasesCounted = 0;
                            
                            for (int i = 0; i < temp.length(); i++) {
                                if (!pDataArray->keepN) {
                                    //eliminate N's
                                    if (toupper(temp[i]) == 'N') {
                                        temp[i] = '.';
                                        tempLength--;
                                        if (tempLength < pDataArray->numbases) { stopSpot = 0; break; }
                                    }
                                }
                                if(isalpha(temp[i])) { numBasesCounted++; }
                                
                                if (numBasesCounted >= pDataArray->numbases) { stopSpot = i; break; }
                            }
                            
                            if (stopSpot == 0) { temp = ""; }
                            else {  temp = temp.substr(0, stopSpot+1);  }
							
                            qualValues = seq.getName() +'\t' + toString(0) + '\t' + toString(stopSpot+1) + '\n';
                            
                        }else { 
                            if (!Short) { temp = ""; qualValues = seq.getName() +'\t' + toString(0) + '\t' + toString(0) + '\n'; } //sequence too short
                            else { qualValues = seq.getName() +'\t' + toString(0) + '\t' + toString(tempLength) + '\n'; }
                        }				
                    }else { //you are keeping the back
                        int tempLength = tempUnaligned.length();
                        if (tempLength > pDataArray->numbases) { //you have enough bases to remove some
                            
                            int stopSpot = 0;
                            int numBasesCounted = 0;
                            
                            for (int i = (temp.length()-1); i >= 0; i--) {
                                if (!pDataArray->keepN) {
                                    //eliminate N's
                                    if (toupper(temp[i]) == 'N') {
                                        temp[i] = '.';
                                        tempLength--;
                                        if (tempLength < pDataArray->numbases) { stopSpot = 0; break; }
                                    }
                                }
                                if(isalpha(temp[i])) { numBasesCounted++; }
                                
                                if (numBasesCounted >= pDataArray->numbases) { stopSpot = i; break; }
                            }
                            
                            if (stopSpot == 0) { temp = ""; }
                            else {  temp = temp.substr(stopSpot);  }
                            
                            qualValues = seq.getName() +'\t' + toString(stopSpot) + '\t' + toString(temp.length()-1) + '\n';
                            
                        }else {
                            if (!Short) { temp = ""; qualValues = seq.getName() +'\t' + toString(0) + '\t' + toString(0) + '\n'; } //sequence too short
                            else { qualValues = seq.getName() +'\t' + toString(0) + '\t' + toString(tempLength) + '\n'; }
                        }
                    }
                }
                
                string newSeqString = temp;
                ///////////////////////////////////////////////////////////////////////
				
				//output trimmed sequence
				if (newSeqString != "") {
					out << ">" << seq.getName() << endl << newSeqString << endl;
				}else{
					outAcc << seq.getName() << endl;
					pDataArray->wroteAccnos = true;
				}
                if (pDataArray->fastaFileTemp != "") {  outfTemp << qualValues << endl;  }
                pDataArray->count++;
			}
            //report progress
			if((pDataArray->count) % 1000 == 0){	pDataArray->m->mothurOut(toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
			
		}
		//report progress
		if((pDataArray->count) % 1000 != 0){	pDataArray->m->mothurOut(toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
        
		
		in.close();
        out.close();
        outAcc.close();
        if (pDataArray->fastaFileTemp != "") { outfTemp.close(); }
				
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ChopsSeqsCommand", "MyChopThreadFunction");
		exit(1);
	}
} 
#endif



#endif


