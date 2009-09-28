/*
 *  chimeraslayer.cpp
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Pschloss Lab. All rights reserved.
 *
 */

#include "chimeraslayer.h"

//***************************************************************************************************************
ChimeraSlayer::ChimeraSlayer(string filename, string temp) {  fastafile = filename;  templateFile = temp;  }
//***************************************************************************************************************

ChimeraSlayer::~ChimeraSlayer() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)			{  delete querySeqs[i];			}
		for (int i = 0; i < templateSeqs.size(); i++)		{  delete templateSeqs[i];		}
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "~ChimeraSlayer");
		exit(1);
	}
}	
//***************************************************************************************************************
void ChimeraSlayer::print(ostream& out) {
	try {
		
		mothurOutEndLine();
		
				
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}

//***************************************************************************************************************
void ChimeraSlayer::getChimeras() {
	try {
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		templateSeqs = readSeqs(templateFile);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		//break up file if needed
		int linesPerProcess = numSeqs / processors ;
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			//find breakup of sequences for all times we will Parallelize
			if (processors == 1) {   lines.push_back(new linePair(0, numSeqs));  }
			else {
				//fill line pairs
				for (int i = 0; i < (processors-1); i++) {			
					lines.push_back(new linePair((i*linesPerProcess), ((i*linesPerProcess) + linesPerProcess)));
				}
				//this is necessary to get remainder of processors / numSeqs so you don't miss any lines at the end
				int i = processors - 1;
				lines.push_back(new linePair((i*linesPerProcess), numSeqs));
			}
		#else
			lines.push_back(new linePair(0, numSeqs));
		#endif
		
		if (seqMask != "") {	decalc = new DeCalculator();	} //to use below
		
		//referenceSeqs, numWanted, matchScore, misMatchPenalty, divR, minSimilarity
		maligner = new Maligner(templateSeqs, numWanted, match, misMatch, 1.01, minSim);
		slayer = new Slayer(window, increment, minSim, divR);
		
		for (int i = 0; i < querySeqs.size(); i++) {
		
			string chimeraFlag = maligner->getResults(querySeqs[i]);
			float percentIdentical = maligner->getPercentID();
			vector<results> Results = maligner->getOutput();
			
			cout << querySeqs[4]->getName() << '\t' << chimeraFlag << '\t' << percentIdentical << endl;
			
			for (int j = 0; j < Results.size(); j++) {
				cout << "regionStart = " << Results[j].regionStart << "\tRegionEnd = " << Results[j].regionEnd << "\tName = " << Results[j].parent << "\tPerQP = " << Results[j].queryToParent << "\tLocalPerQP = " << Results[j].queryToParentLocal << "\tdivR = " << Results[j].divR << endl;
				
			}
			
			if (chimeraFlag == "yes") {
			
				//get sequence that were given from maligner results
				vector<SeqDist> seqs;
				for (int j = 0; j < Results.size(); j++) {
					Sequence* seq = getSequence(Results[j].parent); //makes copy so you can filter and mask and not effect template
					
					//seq = NULL if error occurred in getSequence
					if (seq == NULL) {  break;	}
					else {	
						SeqDist member;
						member.seq = seq;
						member.dist = (Results[j].regionEnd - Results[j].regionStart + 1) * Results[j].queryToParentLocal;
						seqs.push_back(member);	
					}
				}
			
				//limit number of parents to explore - default 5
				if (Results.size() > parents) {
					//sort by distance
					sort(seqs.begin(), seqs.end(), compareSeqDist);
					//prioritize larger more similiar sequence fragments
					reverse(seqs.begin(), seqs.end());
					
					for (int k = seqs.size()-1; k > (parents-1); k--)  {  
						delete seqs[k].seq;
						seqs.pop_back();	
					}
				}
				
				//put seqs into vector to send to slayer
				vector<Sequence*> seqsForSlayer;
				for (int k = 0; k < seqs.size(); k++) {  seqsForSlayer.push_back(seqs[k].seq);	}
			
				//mask then send to slayer...
				if (seqMask != "") {
					decalc->setMask(seqMask);

					//mask querys
					decalc->runMask(querySeqs[i]);
					
					//mask parents
					for (int k = 0; k < seqsForSlayer.size(); k++) {
						decalc->runMask(seqsForSlayer[k]);
					}
					
				}
				
				//send to slayer
				slayer->getResults(querySeqs[i], seqsForSlayer);
				
				//free memory
				for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }
			}
			
		}	
		//free memory
		for (int i = 0; i < lines.size(); i++)					{	delete lines[i];	}
		
		if (seqMask != "") {
			delete decalc; 
		}

			
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence* ChimeraSlayer::getSequence(string name) {
	try{
		Sequence* temp;
		
		//look through templateSeqs til you find it
		int spot = -1;
		for (int i = 0; i < templateSeqs.size(); i++) {
			if (name == templateSeqs[i]->getName()) {  
				spot = i;
				break;
			}
		}
		
		if(spot == -1) { mothurOut("Error: Could not find sequence in chimeraSlayer."); mothurOutEndLine(); return NULL; }
		
		temp = new Sequence(templateSeqs[spot]->getName(), templateSeqs[spot]->getAligned());
		
		return temp;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "getSequence");
		exit(1);
	}
}
//***************************************************************************************************************



