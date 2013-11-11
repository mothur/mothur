#ifndef CHIMERAUCHIMECOMMAND_H
#define CHIMERAUCHIMECOMMAND_H


/*
 *  chimerauchimecommand.h
 *  Mothur
 *
 *  Created by westcott on 5/13/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequenceparser.h"
#include "counttable.h"
#include "sequencecountparser.h"

/***********************************************************/

class ChimeraUchimeCommand : public Command {
public:
	ChimeraUchimeCommand(string);
	ChimeraUchimeCommand();
	~ChimeraUchimeCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.uchime";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "uchime by Robert C. Edgar\nhttp://drive5.com/uchime\nThis code was donated to the public domain.\nEdgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011), UCHIME improves sensitivity and speed of chimera detection.  Bioinformatics 27:2194.\nhttp://www.mothur.org/wiki/Chimera.uchime\n"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	struct linePair {
		int start;
		int end;
		linePair(int i, int j) : start(i), end(j) {}
	};
	
	vector<int> processIDS;   //processid
	int driver(string, string, string, string, int&);
	int createProcesses(string, string, string, string, int&);
		
	bool abort, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, useXa, useChunks, useMinchunk, useIdsmoothwindow, useMinsmoothid, useMaxp, skipgaps, skipgaps2, useMinlen, useMaxlen, ucl, useQueryfract, hasCount, hasName, dups;
	string fastafile, groupfile, templatefile, outputDir, namefile, countfile, abskew, minh, mindiv, xn, dn, xa, chunks, minchunk, idsmoothwindow, minsmoothid, maxp, minlen, maxlen, queryfract, uchimeLocation, strand;
	int processors;
	
	SequenceParser* sparser;
    SequenceCountParser* cparser;
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	vector<string> groupFileNames;
	
	string getNamesFile(string&);
	int readFasta(string, map<string, string>&);
	int printFile(vector<seqPriorityNode>&, string);
	int deconvoluteResults(map<string, string>&, string, string, string);
	int driverGroups(string, string, string, string, string, int, int, vector<string>);
	int createProcessesGroups(string, string, string, string, string, vector<string>, string, string, string);
    int prepFile(string filename, string);


};

/***********************************************************/
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct uchimeData {
	string fastafile; 
	string namefile; 
	string groupfile;
	string outputFName;
	string accnos, alns, filename, templatefile, uchimeLocation, countlist;
	MothurOut* m;
	int start;
	int end;
	int threadID, count, numChimeras;
	vector<string> groups;
	bool dups, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, useXa, useChunks, useMinchunk, useIdsmoothwindow, useMinsmoothid, useMaxp, skipgaps, skipgaps2, useMinlen, useMaxlen, ucl, useQueryfract, hasCount;
	string abskew, minh, mindiv, xn, dn, xa, chunks, minchunk, idsmoothwindow, minsmoothid, maxp, minlen, maxlen, queryfract, strand;
	
	uchimeData(){}
	uchimeData(string o, string uloc, string t, string file, string f, string n, string g, string ac,  string al, string nc, vector<string> gr, MothurOut* mout, int st, int en, int tid) {
		fastafile = f;
		namefile = n;
		groupfile = g;
		filename = file;
		outputFName = o;
		templatefile = t;
		accnos = ac;
		alns = al;
		m = mout;
		start = st;
		end = en;
		threadID = tid;
		groups = gr;
		count = 0;
		numChimeras = 0;
        uchimeLocation = uloc;
        countlist = nc;
	}
	void setBooleans(bool dps, bool Abskew, bool calns, bool MinH, bool Mindiv, bool Xn, bool Dn, bool Xa, bool Chunks, bool Minchunk, bool Idsmoothwindow, bool Minsmoothid, bool Maxp, bool skipgap, bool skipgap2, bool Minlen, bool Maxlen, bool uc, bool Queryfract, bool hc) {
		useAbskew = Abskew;
		chimealns = calns;
		useMinH = MinH;
		useMindiv = Mindiv;
		useXn = Xn;
		useDn = Dn;
		useXa = Xa;
		useChunks = Chunks;
		useMinchunk = Minchunk;
		useIdsmoothwindow = Idsmoothwindow;
		useMinsmoothid = Minsmoothid;
		useMaxp = Maxp;
		skipgaps = skipgap;
		skipgaps2 = skipgap2;
		useMinlen = Minlen;
		useMaxlen = Maxlen;
		ucl = uc;
		useQueryfract = Queryfract;
        hasCount = hc;
        dups = dps;
	}
	
	void setVariables(string abske, string min, string mindi, string x, string d, string xa2, string chunk, string minchun, string idsmoothwindo, string minsmoothi, string max, string minle, string maxle, string queryfrac, string stra) {
		abskew = abske;
		minh = min;
		mindiv = mindi;
        strand = stra;
		xn = x;
		dn = d;
		xa = xa2;
		chunks = chunk;
		minchunk = minchun;
		idsmoothwindow = idsmoothwindo;
		minsmoothid = minsmoothi;
		maxp = max;
		minlen = minle;
		maxlen = maxle;
		queryfract = queryfrac;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyUchimeThreadFunction(LPVOID lpParam){ 
	uchimeData* pDataArray;
	pDataArray = (uchimeData*)lpParam;
	
	try {
		
		pDataArray->outputFName = pDataArray->m->getFullPathName(pDataArray->outputFName);
		pDataArray->filename = pDataArray->m->getFullPathName(pDataArray->filename);
		pDataArray->alns = pDataArray->m->getFullPathName(pDataArray->alns);
		
		//clears files
		ofstream out, out1, out2;
		pDataArray->m->openOutputFile(pDataArray->outputFName, out); out.close(); 
		pDataArray->m->openOutputFile(pDataArray->accnos, out1); out1.close();
		if (pDataArray->chimealns) { pDataArray->m->openOutputFile(pDataArray->alns, out2); out2.close(); }
		
		//parse fasta and name file by group
		SequenceParser* parser;
        SequenceCountParser* cparser;
		if (pDataArray->hasCount) {
            CountTable* ct = new CountTable();
            ct->readTable(pDataArray->namefile, true, false);
            cparser = new SequenceCountParser(pDataArray->fastafile, *ct);
            delete ct;
        }else {
            if (pDataArray->namefile != "") { parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile, pDataArray->namefile);	}
            else							{ parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile);						}
        }
		
		int totalSeqs = 0;
		int numChimeras = 0;
        
        ofstream outCountList;
        if (pDataArray->hasCount && pDataArray->dups) { pDataArray->m->openOutputFile(pDataArray->countlist, outCountList); }

		
		for (int i = pDataArray->start; i < pDataArray->end; i++) {
			int start = time(NULL);	 if (pDataArray->m->control_pressed) {  if (pDataArray->hasCount) { delete cparser; } { delete parser; } return 0; }
			
            
			int error;
            if (pDataArray->hasCount) { 
                error = cparser->getSeqs(pDataArray->groups[i], pDataArray->filename, true); if ((error == 1) || pDataArray->m->control_pressed) {  delete cparser; return 0; }
            }else {
               error = parser->getSeqs(pDataArray->groups[i], pDataArray->filename, true); if ((error == 1) || pDataArray->m->control_pressed) {  delete parser; return 0; } 
            }
			
			//int numSeqs = driver((outputFName + groups[i]), filename, (accnos+ groups[i]), (alns+ groups[i]), numChimeras);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			//to allow for spaces in the path
			string outputFName = "\"" + pDataArray->outputFName+pDataArray->groups[i] + "\"";
			string filename = "\"" + pDataArray->filename + "\"";
			string alns = "\"" + pDataArray->alns+pDataArray->groups[i] + "\"";
			string accnos = pDataArray->accnos+pDataArray->groups[i];
			
			vector<char*> cPara;
			
            string uchimeCommand = pDataArray->uchimeLocation;
            uchimeCommand = "\"" + uchimeCommand + "\"";
			
			char* tempUchime;
			tempUchime= new char[uchimeCommand.length()+1]; 
			*tempUchime = '\0';
			strncat(tempUchime, uchimeCommand.c_str(), uchimeCommand.length());
			cPara.push_back(tempUchime);
			
			char* tempIn = new char[8]; 
			*tempIn = '\0'; strncat(tempIn, "--input", 7);
			//strcpy(tempIn, "--input"); 
			cPara.push_back(tempIn);
			char* temp = new char[filename.length()+1];
			*temp = '\0'; strncat(temp, filename.c_str(), filename.length());
			//strcpy(temp, filename.c_str());
			cPara.push_back(temp);
			
			char* tempO = new char[12]; 
			*tempO = '\0'; strncat(tempO, "--uchimeout", 11);
			//strcpy(tempO, "--uchimeout"); 
			cPara.push_back(tempO);
			char* tempout = new char[outputFName.length()+1];
			//strcpy(tempout, outputFName.c_str());
			*tempout = '\0'; strncat(tempout, outputFName.c_str(), outputFName.length());
			cPara.push_back(tempout);
			
			if (pDataArray->chimealns) {
				char* tempA = new char[13]; 
				*tempA = '\0'; strncat(tempA, "--uchimealns", 12);
				//strcpy(tempA, "--uchimealns"); 
				cPara.push_back(tempA);
				char* tempa = new char[alns.length()+1];
				//strcpy(tempa, alns.c_str());
				*tempa = '\0'; strncat(tempa, alns.c_str(), alns.length());
				cPara.push_back(tempa);
			}
			
            if (pDataArray->strand != "") {
                char* tempA = new char[9]; 
                *tempA = '\0'; strncat(tempA, "--strand", 8);
                cPara.push_back(tempA);
                char* tempa = new char[pDataArray->strand.length()+1];
                *tempa = '\0'; strncat(tempa, pDataArray->strand.c_str(), pDataArray->strand.length());
                cPara.push_back(tempa);
            }
            
			if (pDataArray->useAbskew) {
				char* tempskew = new char[9];
				*tempskew = '\0'; strncat(tempskew, "--abskew", 8);
				//strcpy(tempskew, "--abskew"); 
				cPara.push_back(tempskew);
				char* tempSkew = new char[pDataArray->abskew.length()+1];
				//strcpy(tempSkew, abskew.c_str());
				*tempSkew = '\0'; strncat(tempSkew, pDataArray->abskew.c_str(), pDataArray->abskew.length());
				cPara.push_back(tempSkew);
			}
			
			if (pDataArray->useMinH) {
				char* tempminh = new char[7]; 
				*tempminh = '\0'; strncat(tempminh, "--minh", 6);
				//strcpy(tempminh, "--minh"); 
				cPara.push_back(tempminh);
				char* tempMinH = new char[pDataArray->minh.length()+1];
				*tempMinH = '\0'; strncat(tempMinH, pDataArray->minh.c_str(), pDataArray->minh.length());
				//strcpy(tempMinH, minh.c_str());
				cPara.push_back(tempMinH);
			}
			
			if (pDataArray->useMindiv) {
				char* tempmindiv = new char[9]; 
				*tempmindiv = '\0'; strncat(tempmindiv, "--mindiv", 8);
				//strcpy(tempmindiv, "--mindiv"); 
				cPara.push_back(tempmindiv);
				char* tempMindiv = new char[pDataArray->mindiv.length()+1];
				*tempMindiv = '\0'; strncat(tempMindiv, pDataArray->mindiv.c_str(), pDataArray->mindiv.length());
				//strcpy(tempMindiv, mindiv.c_str());
				cPara.push_back(tempMindiv);
			}
			
			if (pDataArray->useXn) {
				char* tempxn = new char[5]; 
				//strcpy(tempxn, "--xn"); 
				*tempxn = '\0'; strncat(tempxn, "--xn", 4);
				cPara.push_back(tempxn);
				char* tempXn = new char[pDataArray->xn.length()+1];
				//strcpy(tempXn, xn.c_str());
				*tempXn = '\0'; strncat(tempXn, pDataArray->xn.c_str(), pDataArray->xn.length());
				cPara.push_back(tempXn);
			}
			
			if (pDataArray->useDn) {
				char* tempdn = new char[5]; 
				//strcpy(tempdn, "--dn"); 
				*tempdn = '\0'; strncat(tempdn, "--dn", 4);
				cPara.push_back(tempdn);
				char* tempDn = new char[pDataArray->dn.length()+1];
				*tempDn = '\0'; strncat(tempDn, pDataArray->dn.c_str(), pDataArray->dn.length());
				//strcpy(tempDn, dn.c_str());
				cPara.push_back(tempDn);
			}
			
			if (pDataArray->useXa) {
				char* tempxa = new char[5]; 
				//strcpy(tempxa, "--xa"); 
				*tempxa = '\0'; strncat(tempxa, "--xa", 4);
				cPara.push_back(tempxa);
				char* tempXa = new char[pDataArray->xa.length()+1];
				*tempXa = '\0'; strncat(tempXa, pDataArray->xa.c_str(), pDataArray->xa.length());
				//strcpy(tempXa, xa.c_str());
				cPara.push_back(tempXa);
			}
			
			if (pDataArray->useChunks) {
				char* tempchunks = new char[9]; 
				//strcpy(tempchunks, "--chunks"); 
				*tempchunks = '\0'; strncat(tempchunks, "--chunks", 8);
				cPara.push_back(tempchunks);
				char* tempChunks = new char[pDataArray->chunks.length()+1];
				*tempChunks = '\0'; strncat(tempChunks, pDataArray->chunks.c_str(), pDataArray->chunks.length());
				//strcpy(tempChunks, chunks.c_str());
				cPara.push_back(tempChunks);
			}
			
			if (pDataArray->useMinchunk) {
				char* tempminchunk = new char[11]; 
				//strcpy(tempminchunk, "--minchunk"); 
				*tempminchunk = '\0'; strncat(tempminchunk, "--minchunk", 10);
				cPara.push_back(tempminchunk);
				char* tempMinchunk = new char[pDataArray->minchunk.length()+1];
				*tempMinchunk = '\0'; strncat(tempMinchunk, pDataArray->minchunk.c_str(), pDataArray->minchunk.length());
				//strcpy(tempMinchunk, minchunk.c_str());
				cPara.push_back(tempMinchunk);
			}
			
			if (pDataArray->useIdsmoothwindow) {
				char* tempidsmoothwindow = new char[17]; 
				*tempidsmoothwindow = '\0'; strncat(tempidsmoothwindow, "--idsmoothwindow", 16);
				//strcpy(tempidsmoothwindow, "--idsmoothwindow"); 
				cPara.push_back(tempidsmoothwindow);
				char* tempIdsmoothwindow = new char[pDataArray->idsmoothwindow.length()+1];
				*tempIdsmoothwindow = '\0'; strncat(tempIdsmoothwindow, pDataArray->idsmoothwindow.c_str(), pDataArray->idsmoothwindow.length());
				//strcpy(tempIdsmoothwindow, idsmoothwindow.c_str());
				cPara.push_back(tempIdsmoothwindow);
			}
			
			if (pDataArray->useMaxp) {
				char* tempmaxp = new char[7]; 
				//strcpy(tempmaxp, "--maxp"); 
				*tempmaxp = '\0'; strncat(tempmaxp, "--maxp", 6);
				cPara.push_back(tempmaxp);
				char* tempMaxp = new char[pDataArray->maxp.length()+1];
				*tempMaxp = '\0'; strncat(tempMaxp, pDataArray->maxp.c_str(), pDataArray->maxp.length());
				//strcpy(tempMaxp, maxp.c_str());
				cPara.push_back(tempMaxp);
			}
			
			if (!pDataArray->skipgaps) {
				char* tempskipgaps = new char[13]; 
				//strcpy(tempskipgaps, "--[no]skipgaps");
				*tempskipgaps = '\0'; strncat(tempskipgaps, "--noskipgaps", 12);
				cPara.push_back(tempskipgaps);
			}
			
			if (!pDataArray->skipgaps2) {
				char* tempskipgaps2 = new char[14]; 
				//strcpy(tempskipgaps2, "--[no]skipgaps2"); 
				*tempskipgaps2 = '\0'; strncat(tempskipgaps2, "--noskipgaps2", 13);
				cPara.push_back(tempskipgaps2);
			}
			
			if (pDataArray->useMinlen) {
				char* tempminlen = new char[9]; 
				*tempminlen = '\0'; strncat(tempminlen, "--minlen", 8);
				//strcpy(tempminlen, "--minlen"); 
				cPara.push_back(tempminlen);
				char* tempMinlen = new char[pDataArray->minlen.length()+1];
				//strcpy(tempMinlen, minlen.c_str());
				*tempMinlen = '\0'; strncat(tempMinlen, pDataArray->minlen.c_str(), pDataArray->minlen.length());
				cPara.push_back(tempMinlen);
			}
			
			if (pDataArray->useMaxlen) {
				char* tempmaxlen = new char[9]; 
				//strcpy(tempmaxlen, "--maxlen"); 
				*tempmaxlen = '\0'; strncat(tempmaxlen, "--maxlen", 8);
				cPara.push_back(tempmaxlen);
				char* tempMaxlen = new char[pDataArray->maxlen.length()+1];
				*tempMaxlen = '\0'; strncat(tempMaxlen, pDataArray->maxlen.c_str(), pDataArray->maxlen.length());
				//strcpy(tempMaxlen, maxlen.c_str());
				cPara.push_back(tempMaxlen);
			}
			
			if (pDataArray->ucl) {
				char* tempucl = new char[5]; 
				strcpy(tempucl, "--ucl"); 
				cPara.push_back(tempucl);
			}
			
			if (pDataArray->useQueryfract) {
				char* tempqueryfract = new char[13]; 
				*tempqueryfract = '\0'; strncat(tempqueryfract, "--queryfract", 12);
				//strcpy(tempqueryfract, "--queryfract"); 
				cPara.push_back(tempqueryfract);
				char* tempQueryfract = new char[pDataArray->queryfract.length()+1];
				*tempQueryfract = '\0'; strncat(tempQueryfract, pDataArray->queryfract.c_str(), pDataArray->queryfract.length());
				//strcpy(tempQueryfract, queryfract.c_str());
				cPara.push_back(tempQueryfract);
			}
			
			
			char** uchimeParameters;
			uchimeParameters = new char*[cPara.size()];
			string commandString = "";
			for (int j = 0; j < cPara.size(); j++) {  uchimeParameters[j] = cPara[j];  commandString += toString(cPara[j]) + " "; } 
			//int numArgs = cPara.size();
			
			//uchime_main(numArgs, uchimeParameters); 
			//cout << "commandString = " << commandString << endl;
			commandString = "\"" + commandString + "\"";
            
            if (pDataArray->m->debug) { pDataArray->m->mothurOut("[DEBUG]: uchime command = " + commandString + ".\n"); }
            
			system(commandString.c_str());
			
			//free memory
			for(int j = 0; j < cPara.size(); j++)  {  delete cPara[j];  }
			delete[] uchimeParameters; 
			
			//remove "" from filenames
			outputFName = outputFName.substr(1, outputFName.length()-2);
			filename = filename.substr(1, filename.length()-2);
			alns = alns.substr(1, alns.length()-2);
			
			if (pDataArray->m->control_pressed) { if (pDataArray->hasCount) { delete cparser; } { delete parser; } return 0; }
			
			//create accnos file from uchime results
			ifstream in; 
			pDataArray->m->openInputFile(outputFName, in);
			
			ofstream out;
			pDataArray->m->openOutputFile(accnos, out);
            
			
			int num = 0;
			numChimeras = 0;
            map<string, string> thisnamemap;
            map<string, string>::iterator itN;
            if (pDataArray->dups && !pDataArray->hasCount) { thisnamemap = parser->getNameMap(pDataArray->groups[i]); }
                
			while(!in.eof()) {
				
				if (pDataArray->m->control_pressed) { break; }
				
				string name = "";
				string chimeraFlag = "";
				in >> chimeraFlag >> name;
				
				//fix name 
				name = name.substr(0, name.length()-1); //rip off last /
				name = name.substr(0, name.find_last_of('/'));
				
				for (int j = 0; j < 15; j++) {  in >> chimeraFlag; }
				pDataArray->m->gobble(in);
				
				if (chimeraFlag == "Y") {
                    if (pDataArray->dups) {
                        if (!pDataArray->hasCount) { //output redundant names for each group
                            itN = thisnamemap.find(name);
                            if (itN != thisnamemap.end()) {
                                vector<string> tempNames; pDataArray->m->splitAtComma(itN->second, tempNames);
                                for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                            }else { pDataArray->m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); pDataArray->m->control_pressed = true; }

                        }else {
                            out << name << endl;
                            outCountList << name << '\t' << pDataArray->groups[i] << endl;
                        }
                    }else{  out << name << endl;  }
                    numChimeras++;
                }
				num++;
			}
			in.close();
			out.close();
			
			
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			totalSeqs += num;
			pDataArray->numChimeras += numChimeras;
			
			if (pDataArray->m->control_pressed) { if (pDataArray->hasCount) { delete cparser; } { delete parser; } return 0; }
			
			//remove file made for uchime
			pDataArray->m->mothurRemove(filename);
			
			//append files
			pDataArray->m->appendFiles(outputFName, pDataArray->outputFName); pDataArray->m->mothurRemove(outputFName);
			pDataArray->m->appendFiles(accnos, pDataArray->accnos); pDataArray->m->mothurRemove(accnos);
			if (pDataArray->chimealns) { pDataArray->m->appendFiles(alns, pDataArray->alns); pDataArray->m->mothurRemove(alns); }
			
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(num) + " sequences from group " + pDataArray->groups[i] + ".");	pDataArray->m->mothurOutEndLine();					
			
		}	
		
        if (pDataArray->hasCount && pDataArray->dups) { outCountList.close(); }
		pDataArray->count = totalSeqs;
		if (pDataArray->hasCount) { delete cparser; } { delete parser; }
		return totalSeqs;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ChimeraUchimeCommand", "MyUchimeThreadFunction");
		exit(1);
	}
} 
/**************************************************************************************************/

static DWORD WINAPI MyUchimeSeqsThreadFunction(LPVOID lpParam){ 
	uchimeData* pDataArray;
	pDataArray = (uchimeData*)lpParam;
	
	try {
		
		pDataArray->outputFName = pDataArray->m->getFullPathName(pDataArray->outputFName);
		pDataArray->filename = pDataArray->m->getFullPathName(pDataArray->filename);
		pDataArray->alns = pDataArray->m->getFullPathName(pDataArray->alns);
		
		int totalSeqs = 0;
		int numChimeras = 0;
	
		int start = time(NULL);	 if (pDataArray->m->control_pressed) { return 0; }
			
		//to allow for spaces in the path
		string outputFName = "\"" + pDataArray->outputFName + "\"";
		string filename = "\"" + pDataArray->filename + "\"";
		string alns = "\"" + pDataArray->alns+ "\"";
		string templatefile = "\"" + pDataArray->templatefile + "\"";
		string accnos = pDataArray->accnos;
		
		vector<char*> cPara;
		
		string uchimeCommand = pDataArray->uchimeLocation;
        uchimeCommand = "\"" + uchimeCommand + "\"";
        
        char* tempUchime;
        tempUchime= new char[uchimeCommand.length()+1]; 
        *tempUchime = '\0';
        strncat(tempUchime, uchimeCommand.c_str(), uchimeCommand.length());
        cPara.push_back(tempUchime);
		
        string outputFileName = filename.substr(1, filename.length()-2) + ".uchime_formatted";
        //prepFile(filename.substr(1, filename.length()-2), outputFileName);
        //prepFile(filename, outputFileName);
        /******************************************/
        ifstream in23;
        pDataArray->m->openInputFile((filename.substr(1, filename.length()-2)), in23);
        
        ofstream out23;
        pDataArray->m->openOutputFile(outputFileName, out23);
        
        int fcount = 0;
        while (!in23.eof()) {
            if (pDataArray->m->control_pressed) { break;  }
            
            Sequence seq(in23); pDataArray->m->gobble(in23);
            
            if (seq.getName() != "") { seq.printSequence(out23); fcount++; }
        }
        in23.close();
        out23.close();
        /******************************************/
        
        filename = outputFileName;
        filename = "\"" + filename + "\"";
        
        //add reference file
		char* tempRef = new char[5]; 
		//strcpy(tempRef, "--db"); 
		*tempRef = '\0'; strncat(tempRef, "--db", 4);
		cPara.push_back(tempRef);  
		char* tempR = new char[templatefile.length()+1];
		//strcpy(tempR, templatefile.c_str());
		*tempR = '\0'; strncat(tempR, templatefile.c_str(), templatefile.length());
		cPara.push_back(tempR);
        
		char* tempIn = new char[8]; 
		*tempIn = '\0'; strncat(tempIn, "--input", 7);
		//strcpy(tempIn, "--input"); 
		cPara.push_back(tempIn);
		char* temp = new char[filename.length()+1];
		*temp = '\0'; strncat(temp, filename.c_str(), filename.length());
		//strcpy(temp, filename.c_str());
		cPara.push_back(temp);
		
		char* tempO = new char[12]; 
		*tempO = '\0'; strncat(tempO, "--uchimeout", 11);
		//strcpy(tempO, "--uchimeout"); 
		cPara.push_back(tempO);
		char* tempout = new char[outputFName.length()+1];
		//strcpy(tempout, outputFName.c_str());
		*tempout = '\0'; strncat(tempout, outputFName.c_str(), outputFName.length());
		cPara.push_back(tempout);
		
		if (pDataArray->chimealns) {
			char* tempA = new char[13]; 
			*tempA = '\0'; strncat(tempA, "--uchimealns", 12);
			//strcpy(tempA, "--uchimealns"); 
			cPara.push_back(tempA);
			char* tempa = new char[alns.length()+1];
			//strcpy(tempa, alns.c_str());
			*tempa = '\0'; strncat(tempa, alns.c_str(), alns.length());
			cPara.push_back(tempa);
		}
		
        if (pDataArray->strand != "") {
            char* tempA = new char[9]; 
            *tempA = '\0'; strncat(tempA, "--strand", 8);
            cPara.push_back(tempA);
            char* tempa = new char[pDataArray->strand.length()+1];
            *tempa = '\0'; strncat(tempa, pDataArray->strand.c_str(), pDataArray->strand.length());
            cPara.push_back(tempa);
        }
        
		if (pDataArray->useAbskew) {
			char* tempskew = new char[9];
			*tempskew = '\0'; strncat(tempskew, "--abskew", 8);
			//strcpy(tempskew, "--abskew"); 
			cPara.push_back(tempskew);
			char* tempSkew = new char[pDataArray->abskew.length()+1];
			//strcpy(tempSkew, abskew.c_str());
			*tempSkew = '\0'; strncat(tempSkew, pDataArray->abskew.c_str(), pDataArray->abskew.length());
			cPara.push_back(tempSkew);
		}
		
		if (pDataArray->useMinH) {
			char* tempminh = new char[7]; 
			*tempminh = '\0'; strncat(tempminh, "--minh", 6);
			//strcpy(tempminh, "--minh"); 
			cPara.push_back(tempminh);
			char* tempMinH = new char[pDataArray->minh.length()+1];
			*tempMinH = '\0'; strncat(tempMinH, pDataArray->minh.c_str(), pDataArray->minh.length());
			//strcpy(tempMinH, minh.c_str());
			cPara.push_back(tempMinH);
		}
		
		if (pDataArray->useMindiv) {
			char* tempmindiv = new char[9]; 
			*tempmindiv = '\0'; strncat(tempmindiv, "--mindiv", 8);
			//strcpy(tempmindiv, "--mindiv"); 
			cPara.push_back(tempmindiv);
			char* tempMindiv = new char[pDataArray->mindiv.length()+1];
			*tempMindiv = '\0'; strncat(tempMindiv, pDataArray->mindiv.c_str(), pDataArray->mindiv.length());
			//strcpy(tempMindiv, mindiv.c_str());
			cPara.push_back(tempMindiv);
		}
		
		if (pDataArray->useXn) {
			char* tempxn = new char[5]; 
			//strcpy(tempxn, "--xn"); 
			*tempxn = '\0'; strncat(tempxn, "--xn", 4);
			cPara.push_back(tempxn);
			char* tempXn = new char[pDataArray->xn.length()+1];
			//strcpy(tempXn, xn.c_str());
			*tempXn = '\0'; strncat(tempXn, pDataArray->xn.c_str(), pDataArray->xn.length());
			cPara.push_back(tempXn);
		}
		
		if (pDataArray->useDn) {
			char* tempdn = new char[5]; 
			//strcpy(tempdn, "--dn"); 
			*tempdn = '\0'; strncat(tempdn, "--dn", 4);
			cPara.push_back(tempdn);
			char* tempDn = new char[pDataArray->dn.length()+1];
			*tempDn = '\0'; strncat(tempDn, pDataArray->dn.c_str(), pDataArray->dn.length());
			//strcpy(tempDn, dn.c_str());
			cPara.push_back(tempDn);
		}
		
		if (pDataArray->useXa) {
			char* tempxa = new char[5]; 
			//strcpy(tempxa, "--xa"); 
			*tempxa = '\0'; strncat(tempxa, "--xa", 4);
			cPara.push_back(tempxa);
			char* tempXa = new char[pDataArray->xa.length()+1];
			*tempXa = '\0'; strncat(tempXa, pDataArray->xa.c_str(), pDataArray->xa.length());
			//strcpy(tempXa, xa.c_str());
			cPara.push_back(tempXa);
		}
		
		if (pDataArray->useChunks) {
			char* tempchunks = new char[9]; 
			//strcpy(tempchunks, "--chunks"); 
			*tempchunks = '\0'; strncat(tempchunks, "--chunks", 8);
			cPara.push_back(tempchunks);
			char* tempChunks = new char[pDataArray->chunks.length()+1];
			*tempChunks = '\0'; strncat(tempChunks, pDataArray->chunks.c_str(), pDataArray->chunks.length());
			//strcpy(tempChunks, chunks.c_str());
			cPara.push_back(tempChunks);
		}
		
		if (pDataArray->useMinchunk) {
			char* tempminchunk = new char[11]; 
			//strcpy(tempminchunk, "--minchunk"); 
			*tempminchunk = '\0'; strncat(tempminchunk, "--minchunk", 10);
			cPara.push_back(tempminchunk);
			char* tempMinchunk = new char[pDataArray->minchunk.length()+1];
			*tempMinchunk = '\0'; strncat(tempMinchunk, pDataArray->minchunk.c_str(), pDataArray->minchunk.length());
			//strcpy(tempMinchunk, minchunk.c_str());
			cPara.push_back(tempMinchunk);
		}
		
		if (pDataArray->useIdsmoothwindow) {
			char* tempidsmoothwindow = new char[17]; 
			*tempidsmoothwindow = '\0'; strncat(tempidsmoothwindow, "--idsmoothwindow", 16);
			//strcpy(tempidsmoothwindow, "--idsmoothwindow"); 
			cPara.push_back(tempidsmoothwindow);
			char* tempIdsmoothwindow = new char[pDataArray->idsmoothwindow.length()+1];
			*tempIdsmoothwindow = '\0'; strncat(tempIdsmoothwindow, pDataArray->idsmoothwindow.c_str(), pDataArray->idsmoothwindow.length());
			//strcpy(tempIdsmoothwindow, idsmoothwindow.c_str());
			cPara.push_back(tempIdsmoothwindow);
		}
		
		if (pDataArray->useMaxp) {
			char* tempmaxp = new char[7]; 
			//strcpy(tempmaxp, "--maxp"); 
			*tempmaxp = '\0'; strncat(tempmaxp, "--maxp", 6);
			cPara.push_back(tempmaxp);
			char* tempMaxp = new char[pDataArray->maxp.length()+1];
			*tempMaxp = '\0'; strncat(tempMaxp, pDataArray->maxp.c_str(), pDataArray->maxp.length());
			//strcpy(tempMaxp, maxp.c_str());
			cPara.push_back(tempMaxp);
		}
		
		if (!pDataArray->skipgaps) {
			char* tempskipgaps = new char[13]; 
			//strcpy(tempskipgaps, "--[no]skipgaps");
			*tempskipgaps = '\0'; strncat(tempskipgaps, "--noskipgaps", 12);
			cPara.push_back(tempskipgaps);
		}
		
		if (!pDataArray->skipgaps2) {
			char* tempskipgaps2 = new char[14]; 
			//strcpy(tempskipgaps2, "--[no]skipgaps2"); 
			*tempskipgaps2 = '\0'; strncat(tempskipgaps2, "--noskipgaps2", 13);
			cPara.push_back(tempskipgaps2);
		}
		
		if (pDataArray->useMinlen) {
			char* tempminlen = new char[9]; 
			*tempminlen = '\0'; strncat(tempminlen, "--minlen", 8);
			//strcpy(tempminlen, "--minlen"); 
			cPara.push_back(tempminlen);
			char* tempMinlen = new char[pDataArray->minlen.length()+1];
			//strcpy(tempMinlen, minlen.c_str());
			*tempMinlen = '\0'; strncat(tempMinlen, pDataArray->minlen.c_str(), pDataArray->minlen.length());
			cPara.push_back(tempMinlen);
		}
		
		if (pDataArray->useMaxlen) {
			char* tempmaxlen = new char[9]; 
			//strcpy(tempmaxlen, "--maxlen"); 
			*tempmaxlen = '\0'; strncat(tempmaxlen, "--maxlen", 8);
			cPara.push_back(tempmaxlen);
			char* tempMaxlen = new char[pDataArray->maxlen.length()+1];
			*tempMaxlen = '\0'; strncat(tempMaxlen, pDataArray->maxlen.c_str(), pDataArray->maxlen.length());
			//strcpy(tempMaxlen, maxlen.c_str());
			cPara.push_back(tempMaxlen);
		}
		
		if (pDataArray->ucl) {
			char* tempucl = new char[5]; 
			strcpy(tempucl, "--ucl"); 
			cPara.push_back(tempucl);
		}
		
		if (pDataArray->useQueryfract) {
			char* tempqueryfract = new char[13]; 
			*tempqueryfract = '\0'; strncat(tempqueryfract, "--queryfract", 12);
			//strcpy(tempqueryfract, "--queryfract"); 
			cPara.push_back(tempqueryfract);
			char* tempQueryfract = new char[pDataArray->queryfract.length()+1];
			*tempQueryfract = '\0'; strncat(tempQueryfract, pDataArray->queryfract.c_str(), pDataArray->queryfract.length());
			//strcpy(tempQueryfract, queryfract.c_str());
			cPara.push_back(tempQueryfract);
		}
		
		
		char** uchimeParameters;
		uchimeParameters = new char*[cPara.size()];
		string commandString = "";
		for (int j = 0; j < cPara.size(); j++) {  uchimeParameters[j] = cPara[j];  commandString += toString(cPara[j]) + " "; } 
		//int numArgs = cPara.size();
		
        commandString = "\"" + commandString + "\"";
        
		//uchime_main(numArgs, uchimeParameters); 
		//cout << "commandString = " << commandString << endl;
        if (pDataArray->m->debug) { pDataArray->m->mothurOut("[DEBUG]: uchime command = " + commandString + ".\n"); }
		system(commandString.c_str());
		
		//free memory
		for(int j = 0; j < cPara.size(); j++)  {  delete cPara[j];  }
		delete[] uchimeParameters; 
		
		//remove "" from filenames
		outputFName = outputFName.substr(1, outputFName.length()-2);
		filename = filename.substr(1, filename.length()-2);
		alns = alns.substr(1, alns.length()-2);
		
		if (pDataArray->m->control_pressed) { return 0; }
		
		//create accnos file from uchime results
		ifstream in; 
		pDataArray->m->openInputFile(outputFName, in);
		
		ofstream out;
		pDataArray->m->openOutputFile(accnos, out);
		
		numChimeras = 0;
		while(!in.eof()) {
			
			if (pDataArray->m->control_pressed) { break; }
			
			string name = "";
			string chimeraFlag = "";
			in >> chimeraFlag >> name;
			
			for (int j = 0; j < 15; j++) {  in >> chimeraFlag; }
			pDataArray->m->gobble(in);
			
			if (chimeraFlag == "Y") {  out << name << endl; numChimeras++; }
			totalSeqs++;
		}
		in.close();
		out.close();
		
        if (fcount != totalSeqs) { pDataArray->m->mothurOut("[ERROR]: process " + toString(pDataArray->threadID) + " only processed " + toString(pDataArray->count) + " of " + toString(pDataArray->end) + " sequences assigned to it, quitting. \n"); pDataArray->m->control_pressed = true; }
        
		if (pDataArray->m->control_pressed) { return 0; }
		
		pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(totalSeqs) + " sequences.");	pDataArray->m->mothurOutEndLine();					
	
		pDataArray->count = totalSeqs;
		pDataArray->numChimeras = numChimeras;
        
		return totalSeqs;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ChimeraUchimeCommand", "MyUchimeSeqsThreadFunction");
		exit(1);
	}
} 

#endif

/**************************************************************************************************/


#endif


