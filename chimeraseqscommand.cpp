/*
 *  chimeraseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimeraseqscommand.h"
#include "eachgapdist.h"
#include "ignoregaps.h"
#include "onegapdist.h"

//***************************************************************************************************************

ChimeraSeqsCommand::ChimeraSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "filter", "correction", "processors", "method", "window", "increment" };
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; mothurOut("fasta is a required parameter for the chimera.seqs command."); mothurOutEndLine(); abort = true;  }	
			
			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "T"; }
			filter = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "correction", false);		if (temp == "not found") { temp = "T"; }
			correction = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, window);
					
			temp = validParameter.validFile(parameters, "increment", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, increment);
				
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found") { method = "bellerophon"; }
			
			if (method != "bellerophon") { mothurOut(method + " is not a valid method."); mothurOutEndLine();  abort = true; }

		}
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "ChimeraSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraSeqsCommand::help(){
	try {
		mothurOut("The chimera.seqs command reads a fastafile and creates a sorted priority score list of potentially chimeric sequences (ideally, the sequences should already be aligned).\n");
		mothurOut("The chimera.seqs command parameters are fasta, filter, correction, processors and method.  fasta is required.\n");
		mothurOut("The filter parameter allows you to specify if you would like to apply a 50% soft filter.  The default is false. \n");
		mothurOut("The correction parameter allows you to .....  The default is true. \n");
		mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		mothurOut("The method parameter allows you to specify the method for finding chimeric sequences.  The default is bellerophon. \n");
		mothurOut("The chimera.seqs command should be in the following format: \n");
		mothurOut("chimera.seqs(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors, method=bellerophon) \n");
		mothurOut("Example: chimera.seqs(fasta=AD.align, filter=True, correction=true, processors=2, method=yourMethod) \n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "help");
		exit(1);
	}
}
//********************************************************************************************************************
//sorts highest score to lowest
inline bool comparePref(Preference left, Preference right){
	return (left.score[0] > right.score[0]);	
}

//***************************************************************************************************************

ChimeraSeqsCommand::~ChimeraSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		
		//do soft filter
		if (filter)  {
			string optionString = "fasta=" + fastafile + ", soft=50, vertical=F";
			filterSeqs = new FilterSeqsCommand(optionString);
			filterSeqs->execute();
			delete filterSeqs;
			
			//reset fastafile to filtered file
			fastafile = getRootName(fastafile) + "filter.fasta";
		}
		
		distCalculator = new eachGapDist();
		
		//read in sequences
		readSeqs();
		
		int numSeqs = seqs.size();
		
		if (numSeqs == 0) { mothurOut("Error in reading you sequences."); mothurOutEndLine(); return 0; }
		
		//set default window to 25% of sequence length
		string seq0 = seqs[0].getAligned();
		if (window == 0) { window = seq0.length() / 4;  }
		else if (window > (seq0.length() / 2)) {  
			mothurOut("Your sequence length is = " + toString(seq0.length()) + ". You have selected a window size greater than the length of half your aligned sequence. I will run it with a window size of " + toString((seq0.length() / 2))); mothurOutEndLine();
			window = (seq0.length() / 2);
		}
		
		if (increment > (seqs[0].getAlignLength() - (2*window))) { 
			if (increment != 10) {
			
				mothurOut("You have selected a increment that is too large. I will use the default."); mothurOutEndLine();
				increment = 10;
				if (increment > (seqs[0].getAlignLength() - (2*window))) {  increment = 0;  }
				
			}else{ increment = 0; }
		}
cout << "increment = " << increment << endl;		
		if (increment == 0) { iters = 1; }
		else { iters = ((seqs[0].getAlignLength() - (2*window)) / increment); }
		
		//initialize pref
		pref.resize(numSeqs);  
		
		for (int i = 0; i < numSeqs; i++ ) { 
			pref[i].leftParent.resize(2); pref[i].rightParent.resize(2); pref[i].score.resize(2);   pref[i].closestLeft.resize(2); pref[i].closestRight.resize(3);
			pref[i].name = seqs[i].getName();
			pref[i].score[0] = 0.0;  pref[i].score[1] = 0.0; 
			pref[i].closestLeft[0] = 100000.0;  pref[i].closestLeft[1] = 100000.0;  
			pref[i].closestRight[0] = 100000.0;  pref[i].closestRight[1] = 100000.0;  
		}

		int midpoint = window;
		int count = 0;
		while (count < iters) {
				
				//create 2 vectors of sequences, 1 for left side and one for right side
				vector<Sequence> left;  vector<Sequence> right;
				
				for (int i = 0; i < seqs.size(); i++) {
//cout << "whole = " << seqs[i].getAligned() << endl;
					//save left side
					string seqLeft = seqs[i].getAligned().substr(midpoint-window, window);
					Sequence tempLeft;
					tempLeft.setName(seqs[i].getName());
					tempLeft.setAligned(seqLeft);
					left.push_back(tempLeft);
//cout << "left = " << tempLeft.getAligned() << endl;			
					//save right side
					string seqRight = seqs[i].getAligned().substr(midpoint, window);
					Sequence tempRight;
					tempRight.setName(seqs[i].getName());
					tempRight.setAligned(seqRight);
					right.push_back(tempRight);
//cout << "right = " << seqRight << endl;	
				}
				
				//adjust midpoint by increment
				midpoint += increment;
				
				
				//this should be parallelized
				//perference = sum of (| distance of my left to sequence j's left - distance of my right to sequence j's right | )
				//create a matrix containing the distance from left to left and right to right
				//calculate distances
				SparseMatrix* SparseLeft = new SparseMatrix();
				SparseMatrix* SparseRight = new SparseMatrix();
				
				createSparseMatrix(0, left.size(), SparseLeft, left);
				createSparseMatrix(0, right.size(), SparseRight, right);
				
				vector<SeqMap> distMapRight;
				vector<SeqMap> distMapLeft;
				
				// Create a data structure to quickly access the distance information.
				// It consists of a vector of distance maps, where each map contains
				// all distances of a certain sequence. Vector and maps are accessed
				// via the index of a sequence in the distance matrix
				distMapRight = vector<SeqMap>(numSeqs); 
				distMapLeft = vector<SeqMap>(numSeqs); 
				//cout << "left" << endl << endl;
				for (MatData currentCell = SparseLeft->begin(); currentCell != SparseLeft->end(); currentCell++) {
					distMapLeft[currentCell->row][currentCell->column] = currentCell->dist;
					//cout << " i = " << currentCell->row << " j = " << currentCell->column << " dist = " << currentCell->dist << endl;
				}
				//cout << "right" << endl << endl;
				for (MatData currentCell = SparseRight->begin(); currentCell != SparseRight->end(); currentCell++) {
					distMapRight[currentCell->row][currentCell->column] = currentCell->dist;
					//cout << " i = " << currentCell->row << " j = " << currentCell->column << " dist = " << currentCell->dist << endl;
				}
				
				delete SparseLeft;
				delete SparseRight;
				
				
				//fill preference structure
				generatePreferences(distMapLeft, distMapRight, midpoint);
				
				count++;
				
		}
		
		delete distCalculator;
		
		//rank preference score to eachother
		float dme = 0.0;
		float expectedPercent = 1 / (float) (pref.size());
		
		for (int i = 0; i < pref.size(); i++) {	 dme += pref[i].score[0];  }
	
		for (int i = 0; i < pref.size(); i++) {

			//gives the actual percentage of the dme this seq adds
			pref[i].score[0] = pref[i].score[0] / dme;
			
			//how much higher or lower is this than expected
			pref[i].score[0] = pref[i].score[0] / expectedPercent;
			
		}
		
		
		//sort Preferences highest to lowest
		sort(pref.begin(), pref.end(), comparePref);
		
		string outputFileName = getRootName(fastafile) + "chimeras";
		ofstream out;
		openOutputFile(outputFileName, out);
		
		int above1 = 0;
		out << "Name\tScore\tLeft\tRight\t" << endl;
		//output prefenence structure to .chimeras file
		for (int i = 0; i < pref.size(); i++) {
			out << pref[i].name << '\t' << pref[i].score[0] << '\t' << pref[i].leftParent[0] << '\t' << pref[i].rightParent[0] << endl;
			
			//calc # of seqs with preference above 1.0
			if (pref[i].score[0] > 1.0) { 
				above1++; 
				mothurOut(pref[i].name + " is a suspected chimera at breakpoint " + toString(pref[i].midpoint)); mothurOutEndLine();
				mothurOut("It's score is " + toString(pref[i].score[0]) + " with suspected left parent " + pref[i].leftParent[0] + " and right parent " + pref[i].rightParent[0]); mothurOutEndLine();
			}
			
			
		}
		
		//output results to screen
		mothurOutEndLine();
		mothurOut("Sequence with preference score above 1.0: " + toString(above1)); mothurOutEndLine();
		int spot;
		spot = pref.size()-1;
		mothurOut("Minimum:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.975;
		mothurOut("2.5%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.75;
		mothurOut("25%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.50;
		mothurOut("Median: \t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.25;
		mothurOut("75%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.025;
		mothurOut("97.5%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = 0;
		mothurOut("Maximum:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************
void ChimeraSeqsCommand::readSeqs(){
	try {
		ifstream inFASTA;
		openInputFile(fastafile, inFASTA);
		
		//read in seqs and store in vector
		while(!inFASTA.eof()){
			Sequence current(inFASTA);
			
			if (current.getAligned() == "") { current.setAligned(current.getUnaligned()); }
			
			seqs.push_back(current);
			
			gobble(inFASTA);
		}
		inFASTA.close();

	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "readSeqs");
		exit(1);
	}
}

/***************************************************************************************************************/
int ChimeraSeqsCommand::createSparseMatrix(int startSeq, int endSeq, SparseMatrix* sparse, vector<Sequence> s){
	try {

		for(int i=startSeq; i<endSeq; i++){
			
			for(int j=0;j<i;j++){
			
				distCalculator->calcDist(s[i], s[j]);
				float dist = distCalculator->getDist();
			
				PCell temp(i, j, dist);
				sparse->addCell(temp);
				
			}
		}
			
	
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "createSparseMatrix");
		exit(1);
	}
}
/***************************************************************************************************************/
void ChimeraSeqsCommand::generatePreferences(vector<SeqMap> left, vector<SeqMap> right, int mid){
	try {
		
		float dme = 0.0;
		SeqMap::iterator itR;
		SeqMap::iterator itL;
		
		//initialize pref[i]
		for (int i = 0; i < pref.size(); i++) {
			pref[i].score[1] = 0.0;
			pref[i].closestLeft[1] = 100000.0; 
			pref[i].closestRight[1] = 100000.0; 
			pref[i].leftParent[1] = "";
			pref[i].rightParent[1] = "";
		}
	
		for (int i = 0; i < left.size(); i++) {
			
			SeqMap currentLeft = left[i];    //example i = 3;   currentLeft is a map of 0 to the distance of sequence 3 to sequence 0,
												//										1 to the distance of sequence 3 to sequence 1,
												//										2 to the distance of sequence 3 to sequence 2.
			SeqMap currentRight = right[i];		// same as left but with distances on the right side.
			
			for (int j = 0; j < i; j++) {
				
				itL = currentLeft.find(j);
				itR = currentRight.find(j);
//cout << " i = " << i << " j = " << j << " distLeft = " << itL->second << endl;
//cout << " i = " << i << " j = " << j << " distright = " << itR->second << endl;
				
				//if you can find this entry update the preferences
				if ((itL != currentLeft.end()) && (itR != currentRight.end())) {
				
					if (!correction) {
						pref[i].score[1] += abs((itL->second - itR->second));
						pref[j].score[1] += abs((itL->second - itR->second));
//cout << "left " << i << " " << j << " = " << itL->second << " right " << i << " " << j << " = " << itR->second << endl;
//cout << "abs = " << abs((itL->second - itR->second)) << endl;
//cout << i << " score = " << pref[i].score[1] << endl;
//cout << j << " score = " << pref[j].score[1] << endl;
					}else {
						pref[i].score[1] += abs((sqrt(itL->second) - sqrt(itR->second)));
						pref[j].score[1] += abs((sqrt(itL->second) - sqrt(itR->second)));
//cout << "left " << i << " " << j << " = " << itL->second << " right " << i << " " << j << " = " << itR->second << endl;
//cout << "abs = " << abs((sqrt(itL->second) - sqrt(itR->second))) << endl;
//cout << i << " score = " << pref[i].score[1] << endl;
//cout << j << " score = " << pref[j].score[1] << endl;
					}
//cout << "pref[" << i << "].closestLeft[1] = "	<< 	pref[i].closestLeft[1] << " parent = " << pref[i].leftParent[1] << endl;			
					//are you the closest left sequence
					if (itL->second < pref[i].closestLeft[1]) {  

						pref[i].closestLeft[1] = itL->second;
						pref[i].leftParent[1] = seqs[j].getName();
//cout << "updating closest left to " << pref[i].leftParent[1] << endl;
					}
//cout << "pref[" << j << "].closestLeft[1] = "	<< 	pref[j].closestLeft[1] << " parent = " << pref[j].leftParent[1] << endl;	
					if (itL->second < pref[j].closestLeft[1]) { 
						pref[j].closestLeft[1] = itL->second;
						pref[j].leftParent[1] = seqs[i].getName();
//cout << "updating closest left to " << pref[j].leftParent[1] << endl;
					}
					
					//are you the closest right sequence
					if (itR->second < pref[i].closestRight[1]) {   
						pref[i].closestRight[1] = itR->second;
						pref[i].rightParent[1] = seqs[j].getName();
					}
					if (itR->second < pref[j].closestRight[1]) {   
						pref[j].closestRight[1] = itR->second;
						pref[j].rightParent[1] = seqs[i].getName();
					}
					
				}
			}
		
		}
		
		
		  
		//calculate the dme
		
		int count0 = 0;
		for (int i = 0; i < pref.size(); i++) {	 dme += pref[i].score[1];  if (pref[i].score[1] == 0.0) { count0++; }  }
		
		float expectedPercent = 1 / (float) (pref.size() - count0);
//cout << endl << "dme = " << dme << endl;
		//recalculate prefernences based on dme
		for (int i = 0; i < pref.size(); i++) {
//cout << "unadjusted pref " << i << " = " << pref[i].score[1] << endl;	
			// gives the actual percentage of the dme this seq adds
			pref[i].score[1] = pref[i].score[1] / dme;
			
			//how much higher or lower is this than expected
			pref[i].score[1] = pref[i].score[1] / expectedPercent;
			
			//pref[i].score[1] = dme / (dme - 2 * pref[i].score[1]);
			
			//so a non chimeric sequence would be around 1, and a chimeric would be signifigantly higher.
//cout << "adjusted pref " << i << " = " << pref[i].score[1] << endl;					
		}
		
		//is this score bigger then the last score
		for (int i = 0; i < pref.size(); i++) {	 
			
			//update biggest score
			if (pref[i].score[1] > pref[i].score[0]) {
				pref[i].score[0] = pref[i].score[1];
				pref[i].leftParent[0] = pref[i].leftParent[1];
				pref[i].rightParent[0] = pref[i].rightParent[1];
				pref[i].closestLeft[0] = pref[i].closestLeft[1];
				pref[i].closestRight[0] = pref[i].closestRight[1];
				pref[i].midpoint = mid;
			}
			
		}

	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "generatePreferences");
		exit(1);
	}
}
/**************************************************************************************************/

/**************************************************************************************************/

