#ifndef CHIMERAPERSEUSCOMMAND_H
#define CHIMERAPERSEUSCOMMAND_H


/*
 *  chimeraperseuscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/26/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */



#include "mothur.h"
#include "command.hpp"
#include "sequenceparser.h"
#include "myPerseus.h"

/***********************************************************/
class ChimeraPerseusCommand : public Command {
public:
	ChimeraPerseusCommand(string);
	ChimeraPerseusCommand();
	~ChimeraPerseusCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.perseus";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Chimera.perseus\n"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	struct linePair {
		int start;
		int end;
		linePair(int i, int j) : start(i), end(j) {}
	};
	
	bool abort;
	string fastafile, groupfile, outputDir, namefile;
	int processors;
	double cutoff, alpha, beta;
	
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	vector<string> groupFileNames;
	
	string getNamesFile(string&);
	int driver(string, vector<seqData>&, string, int&);
	vector<seqData> readFiles(string, string);
	vector<seqData> loadSequences(SequenceParser&, string);
	int deconvoluteResults(SequenceParser&, string, string);
	int driverGroups(SequenceParser&, string, string, int, int, vector<string>);
	int createProcessesGroups(SequenceParser&, string, string, vector<string>, string, string, string);
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct perseusData {
	string fastafile; 
	string namefile; 
	string groupfile;
	string outputFName;
	string accnos;
	MothurOut* m;
	int start;
	int end;
	int threadID, count, numChimeras;
	double alpha, beta, cutoff;
	vector<string> groups;
	
	perseusData(){}
	perseusData(double a, double b, double c, string o,  string f, string n, string g, string ac, vector<string> gr, MothurOut* mout, int st, int en, int tid) {
		alpha = a;
		beta = b;
		cutoff = c;
		fastafile = f;
		namefile = n;
		groupfile = g;
		outputFName = o;
		accnos = ac;
		m = mout;
		start = st;
		end = en;
		threadID = tid;
		groups = gr;
		count = 0;
		numChimeras = 0;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
#else
static DWORD WINAPI MyPerseusThreadFunction(LPVOID lpParam){ 
	perseusData* pDataArray;
	pDataArray = (perseusData*)lpParam;
	
	try {
		
		//clears files
		ofstream out, out1, out2;
		pDataArray->m->openOutputFile(pDataArray->outputFName, out); out.close(); 
		pDataArray->m->openOutputFile(pDataArray->accnos, out1); out1.close();
		
		//parse fasta and name file by group
		SequenceParser* parser;
		if (pDataArray->namefile != "") { parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile, pDataArray->namefile);	}
		else							{ parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile);						}
		
		int totalSeqs = 0;
		int numChimeras = 0;
		
		for (int i = pDataArray->start; i < pDataArray->end; i++) {
			
			int start = time(NULL);	 if (pDataArray->m->control_pressed) {  delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); return 0; }
			
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("Checking sequences from group " + pDataArray->groups[i] + "...");	pDataArray->m->mothurOutEndLine();					
			
			//vector<seqData> sequences = loadSequences(parser, groups[i]); - same function below
			////////////////////////////////////////////////////////////////////////////////////////
			vector<Sequence> thisGroupsSeqs = parser->getSeqs(pDataArray->groups[i]);
			map<string, string> nameMap = parser->getNameMap(pDataArray->groups[i]);
			map<string, string>::iterator it;
			
			vector<seqData> sequences;
			bool error = false;
			
			for (int j = 0; j < thisGroupsSeqs.size(); j++) {
				
				if (pDataArray->m->control_pressed) {  delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); return 0; }
				
				it = nameMap.find(thisGroupsSeqs[j].getName());
				if (it == nameMap.end()) { error = true; pDataArray->m->mothurOut("[ERROR]: " + thisGroupsSeqs[j].getName() + " is in your fasta file and not in your namefile, please correct."); pDataArray->m->mothurOutEndLine(); }
				else {
					int num = pDataArray->m->getNumNames(it->second);
					sequences.push_back(seqData(thisGroupsSeqs[j].getName(), thisGroupsSeqs[j].getUnaligned(), num));
				}
			}
			
			if (error) { pDataArray->m->control_pressed = true; }
			
			//sort by frequency
			sort(sequences.rbegin(), sequences.rend());
			////////////////////////////////////////////////////////////////////////////////////////

			if (pDataArray->m->control_pressed) { delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); return 0; }
			
			//int numSeqs = driver((outputFName + groups[i]), sequences, (accnos+groups[i]), numChimeras); - same function below
			////////////////////////////////////////////////////////////////////////////////////////
			string chimeraFileName = pDataArray->outputFName+pDataArray->groups[i];
			string accnosFileName = pDataArray->accnos+pDataArray->groups[i];
			
			vector<vector<double> > correctModel(4);	//could be an option in the future to input own model matrix
			for(int j=0;j<4;j++){	correctModel[j].resize(4);	}
			
			correctModel[0][0] = 0.000000;	//AA
			correctModel[1][0] = 11.619259;	//CA
			correctModel[2][0] = 11.694004;	//TA
			correctModel[3][0] = 7.748623;	//GA
			
			correctModel[1][1] = 0.000000;	//CC
			correctModel[2][1] = 7.619657;	//TC
			correctModel[3][1] = 12.852562;	//GC
			
			correctModel[2][2] = 0.000000;	//TT
			correctModel[3][2] = 10.964048;	//TG
			
			correctModel[3][3] = 0.000000;	//GG
			
			for(int k=0;k<4;k++){
				for(int j=0;j<k;j++){
					correctModel[j][k] = correctModel[k][j];
				}
			}
			
			int numSeqs = sequences.size();
			int alignLength = sequences[0].sequence.size();
			
			ofstream chimeraFile;
			ofstream accnosFile;
			pDataArray->m->openOutputFile(chimeraFileName, chimeraFile); 
			pDataArray->m->openOutputFile(accnosFileName, accnosFile); 
			
			Perseus myPerseus;
			vector<vector<double> > binMatrix = myPerseus.binomial(alignLength);
			
			chimeraFile << "SequenceIndex\tName\tDiffsToBestMatch\tBestMatchIndex\tBestMatchName\tDiffstToChimera\tIndexofLeftParent\tIndexOfRightParent\tNameOfLeftParent\tNameOfRightParent\tDistanceToBestMatch\tcIndex\t(cIndex - singleDist)\tloonIndex\tMismatchesToChimera\tMismatchToTrimera\tChimeraBreakPoint\tLogisticProbability\tTypeOfSequence\n";
			
			vector<bool> chimeras(numSeqs, 0);
			
			for(int j=0;j<numSeqs;j++){	
				
				if (pDataArray->m->control_pressed) { delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
				
				vector<bool> restricted = chimeras;
				
				vector<vector<int> > leftDiffs(numSeqs);
				vector<vector<int> > leftMaps(numSeqs);
				vector<vector<int> > rightDiffs(numSeqs);
				vector<vector<int> > rightMaps(numSeqs);
				
				vector<int> singleLeft, bestLeft;
				vector<int> singleRight, bestRight;
				
				int bestSingleIndex, bestSingleDiff;
				vector<pwAlign> alignments(numSeqs);
				
				int comparisons = myPerseus.getAlignments(j, sequences, alignments, leftDiffs, leftMaps, rightDiffs, rightMaps, bestSingleIndex, bestSingleDiff, restricted);
				
				if (pDataArray->m->control_pressed) { delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
				
				int minMismatchToChimera, leftParentBi, rightParentBi, breakPointBi;
				
				string dummyA, dummyB;
				
				if(comparisons >= 2){	
					minMismatchToChimera = myPerseus.getChimera(sequences, leftDiffs, rightDiffs, leftParentBi, rightParentBi, breakPointBi, singleLeft, bestLeft, singleRight, bestRight, restricted);
					
					if (pDataArray->m->control_pressed) { delete parser;  pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
					
					int minMismatchToTrimera = numeric_limits<int>::max();
					int leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB;
					
					if(minMismatchToChimera >= 3 && comparisons >= 3){
						minMismatchToTrimera = myPerseus.getTrimera(sequences, leftDiffs, leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB, singleLeft, bestLeft, singleRight, bestRight, restricted);
						
						if (pDataArray->m->control_pressed) { delete parser;  pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
					}
					
					double singleDist = myPerseus.modeledPairwiseAlignSeqs(sequences[j].sequence, sequences[bestSingleIndex].sequence, dummyA, dummyB, correctModel);
					
					if (pDataArray->m->control_pressed) { delete parser;  pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
					
					string type;
					string chimeraRefSeq;
					
					if(minMismatchToChimera - minMismatchToTrimera >= 3){
						type = "trimera";
						chimeraRefSeq = myPerseus.stitchTrimera(alignments, leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB, leftMaps, rightMaps);
					}
					else{
						type = "chimera";
						chimeraRefSeq = myPerseus.stitchBimera(alignments, leftParentBi, rightParentBi, breakPointBi, leftMaps, rightMaps);
					}
					
					if (pDataArray->m->control_pressed) { delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
					
					double chimeraDist = myPerseus.modeledPairwiseAlignSeqs(sequences[j].sequence, chimeraRefSeq, dummyA, dummyB, correctModel);
					
					if (pDataArray->m->control_pressed) { delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
					
					double cIndex = chimeraDist;//modeledPairwiseAlignSeqs(sequences[j].sequence, chimeraRefSeq);
					double loonIndex = myPerseus.calcLoonIndex(sequences[j].sequence, sequences[leftParentBi].sequence, sequences[rightParentBi].sequence, breakPointBi, binMatrix);		
					
					if (pDataArray->m->control_pressed) { delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); chimeraFile.close(); pDataArray->m->mothurRemove(chimeraFileName); accnosFile.close(); pDataArray->m->mothurRemove(accnosFileName); return 0; }
					
					chimeraFile << j << '\t' << sequences[j].seqName << '\t' << bestSingleDiff << '\t' << bestSingleIndex << '\t' << sequences[bestSingleIndex].seqName << '\t';
					chimeraFile << minMismatchToChimera << '\t' << leftParentBi << '\t' << rightParentBi << '\t' << sequences[leftParentBi].seqName << '\t' << sequences[rightParentBi].seqName << '\t';
					chimeraFile << singleDist << '\t' << cIndex << '\t' << (cIndex - singleDist) << '\t' << loonIndex << '\t';
					chimeraFile << minMismatchToChimera << '\t' << minMismatchToTrimera << '\t' << breakPointBi << '\t';
					
					double probability = myPerseus.classifyChimera(singleDist, cIndex, loonIndex, pDataArray->alpha, pDataArray->beta);
					
					chimeraFile << probability << '\t';
					
					if(probability > pDataArray->cutoff){ 
						chimeraFile << type << endl;
						accnosFile << sequences[j].seqName << endl;
						chimeras[j] = 1;
						numChimeras++;
					}
					else{
						chimeraFile << "good" << endl;
					}
					
				}
				else{
					chimeraFile << j << '\t' << sequences[j].seqName << "\t0\t0\tNull\t0\t0\t0\tNull\tNull\t0.0\t0.0\t0.0\t0\t0\t0\t0.0\t0.0\tgood" << endl;
				}
				//report progress
				if((j+1) % 100 == 0){ 	pDataArray->m->mothurOut("Processing sequence: " + toString(j+1) + "\n");		}
			}
			
			if((numSeqs) % 100 != 0){ 	pDataArray->m->mothurOut("Processing sequence: " + toString(numSeqs) + "\n");		}
			
			chimeraFile.close();
			accnosFile.close();
			////////////////////////////////////////////////////////////////////////////////////////

			totalSeqs += numSeqs;
			
			//append files
			pDataArray->m->appendFiles(chimeraFileName, pDataArray->outputFName); pDataArray->m->mothurRemove(chimeraFileName);
			pDataArray->m->appendFiles(accnosFileName, pDataArray->accnos); pDataArray->m->mothurRemove(accnosFileName);
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences from group " + pDataArray->groups[i] + ".");	pDataArray->m->mothurOutEndLine();					
			
			if (pDataArray->m->control_pressed) { delete parser; pDataArray->m->mothurRemove(pDataArray->outputFName); pDataArray->m->mothurRemove(pDataArray->accnos); return 0; }
		}	
		
		pDataArray->count = totalSeqs;
		delete parser;
		return totalSeqs;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ChimeraUchimeCommand", "MyPerseusThreadFunction");
		exit(1);
	}
} 
/**************************************************************************************************/

#endif

#endif


