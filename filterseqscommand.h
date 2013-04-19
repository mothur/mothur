#ifndef FILTERSEQSCOMMAND_H
#define FILTERSEQSCOMMAND_H

/*
 *  filterseqscommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "filters.h"

class Sequence;
class FilterSeqsCommand : public Command {

public:
	FilterSeqsCommand(string);
	FilterSeqsCommand();
	~FilterSeqsCommand() {};
	
	vector<string> setParameters();
	string getCommandName()			{ return "filter.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Filter.seqs"; }
	string getDescription()		{ return "removes columns from alignments based on a criteria defined by the user"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};

	vector<linePair*> lines;
	vector<int> processIDS;
    map<int, vector<unsigned long long> > savedPositions;

	string vertical, filter, fasta, hard, outputDir, filterFileName;
	vector<string> fastafileNames;	
	int alignmentLength, processors;
	vector<int> bufferSizes;
	vector<string> outputNames;

	char trump;
	bool abort;
	float soft;
	int numSeqs;
	
	string createFilter();
	int filterSequences();
	int createProcessesCreateFilter(Filters&, string);
	int createProcessesRunFilter(string, string, string);
	int driverRunFilter(string, string, string, linePair*);
	int driverCreateFilter(Filters& F, string filename, linePair* line);
	#ifdef USE_MPI
	int driverMPIRun(int, int, MPI_File&, MPI_File&, vector<unsigned long long>&);
	int MPICreateFilter(int, int, Filters&, MPI_File&, vector<unsigned long long>&);	
	#endif
	
};


/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct filterData {
	Filters F;
    int count, tid, alignmentLength;
    unsigned long long start, end;
    MothurOut* m;
    string filename, vertical, hard;
    char trump;
    float soft;
	
	filterData(){}
	filterData(string fn, MothurOut* mout, unsigned long long st, unsigned long long en, int aLength, char tr, string vert, float so, string ha, int t) {
        filename = fn;
		m = mout;
		start = st;
		end = en;
        tid = t;
        trump = tr;
        alignmentLength = aLength;
        vertical = vert;
        soft = so;
        hard = ha;
		count = 0;
	}
};
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct filterRunData {
    int count, tid, alignmentLength;
    unsigned long long start, end;
    MothurOut* m;
    string filename;
    string filter, outputFilename;
	
	filterRunData(){}
	filterRunData(string f, string fn, string ofn, MothurOut* mout, unsigned long long st, unsigned long long en, int aLength, int t) {
        filter = f;
        outputFilename = ofn;
        filename = fn;
		m = mout;
		start = st;
		end = en;
        tid = t;
        alignmentLength = aLength;
		count = 0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyCreateFilterThreadFunction(LPVOID lpParam){ 
	filterData* pDataArray;
	pDataArray = (filterData*)lpParam;
	
	try {

		if (pDataArray->soft != 0)			{  pDataArray->F.setSoft(pDataArray->soft);		}
		if (pDataArray->trump != '*')		{  pDataArray->F.setTrump(pDataArray->trump);	}
		
		pDataArray->F.setLength(pDataArray->alignmentLength);
		
		if(pDataArray->trump != '*' || pDataArray->m->isTrue(pDataArray->vertical) || pDataArray->soft != 0){
			pDataArray->F.initialize();
		}
		
		if(pDataArray->hard.compare("") != 0)	{	pDataArray->F.doHard(pDataArray->hard);		}
		else						{	pDataArray->F.setFilter(string(pDataArray->alignmentLength, '1'));	}
        
		ifstream in;
		pDataArray->m->openInputFile(pDataArray->filename, in);
        
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			in.seekg(0);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}
		
		pDataArray->count = 0;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) { in.close(); pDataArray->count = 1; return 1; }
			
			Sequence current(in); pDataArray->m->gobble(in); 
			
			if (current.getName() != "") {
				if (current.getAligned().length() != pDataArray->alignmentLength) { pDataArray->m->mothurOut("Sequences are not all the same length, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true;  }
                
                if(pDataArray->trump != '*')			{	pDataArray->F.doTrump(current);		}
                if(pDataArray->m->isTrue(pDataArray->vertical) || pDataArray->soft != 0)	{	pDataArray->F.getFreqs(current);	}
			}
            pDataArray->count++;
            //report progress
			if((i) % 100 == 0){	pDataArray->m->mothurOutJustToScreen(toString(i)+"\n"); 		}
		}
		
        if((pDataArray->count) % 100 != 0){	pDataArray->m->mothurOutJustToScreen(toString(pDataArray->count)+"\n"); 		}
        
		in.close();
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "FilterSeqsCommand", "MyCreateFilterThreadFunction");
		exit(1);
	}
} 
/**************************************************************************************************/
static DWORD WINAPI MyRunFilterThreadFunction(LPVOID lpParam){ 
	filterRunData* pDataArray;
	pDataArray = (filterRunData*)lpParam;
	
	try {
        
        ofstream out;
		pDataArray->m->openOutputFile(pDataArray->outputFilename, out);

		ifstream in;
		pDataArray->m->openInputFile(pDataArray->filename, in);
        
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			in.seekg(0);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}
		
		pDataArray->count = 0;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) { in.close(); out.close(); pDataArray->count = 1; return 1; }
			
			Sequence seq(in); pDataArray->m->gobble(in);
            if (seq.getName() != "") {
                string align = seq.getAligned();
                string filterSeq = "";
                
                for(int j=0;j<pDataArray->alignmentLength;j++){
                    if(pDataArray->filter[j] == '1'){
                        filterSeq += align[j];
                    }
                }
                
                out << '>' << seq.getName() << endl << filterSeq << endl;
            }
            pDataArray->count++;
            //report progress
			if((i) % 100 == 0){	pDataArray->m->mothurOutJustToScreen(toString(i)+"\n"); 		}
		}
		
        if((pDataArray->count) % 100 != 0){	pDataArray->m->mothurOutJustToScreen(toString(pDataArray->count)+"\n"); 		}
        
		in.close();
        out.close();
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "FilterSeqsCommand", "MyRunFilterThreadFunction");
		exit(1);
	}
} 
/**************************************************************************************************/
#endif


#endif
