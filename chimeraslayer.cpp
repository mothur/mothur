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
		mothurOut("Only reporting sequence supported by 90% of bootstrapped results.");
		mothurOutEndLine();
		
		for (int i = 0; i < querySeqs.size(); i++) {
		
			if (chimeraFlags[i] == "yes") {	
				if ((chimeraResults[i][0].bsa >= 90.0) || (chimeraResults[i][0].bsb >= 90.0)) {
					mothurOut(querySeqs[i]->getName() + "\tyes"); mothurOutEndLine();
					out << querySeqs[i]->getName() << "\tyes" << endl;
				}else {
					out << querySeqs[i]->getName() << "\tno" << endl;
					mothurOut(querySeqs[i]->getName() + "\tno"); mothurOutEndLine();
				}

				printBlock(chimeraResults[i][0], out, i);
				
				out << endl;
			}else{
				out << querySeqs[i]->getName() << "\tno" << endl;
				mothurOut(querySeqs[i]->getName() + "\tno"); mothurOutEndLine();
			}
		}
				
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
		
		chimeraResults.resize(numSeqs);
		chimeraFlags.resize(numSeqs, "no");
		spotMap.resize(numSeqs);
		
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
		
		//initialize spotMap
		for (int j = 0; j < numSeqs; j++) {
			for (int i = 0; i < querySeqs[0]->getAligned().length(); i++) {
				spotMap[j][i] = i;
			}
		}
		
		//referenceSeqs, numWanted, matchScore, misMatchPenalty, divR, minSimilarity
		maligner = new Maligner(templateSeqs, numWanted, match, misMatch, 1.01, minSim);
		slayer = new Slayer(window, increment, minSim, divR, iters);
		
		for (int i = 0; i < querySeqs.size(); i++) {
		
			string chimeraFlag = maligner->getResults(querySeqs[i]);
			//float percentIdentical = maligner->getPercentID();
			vector<results> Results = maligner->getOutput();
			
			cout << "Processing sequence: " << i+1 << endl;
			
			//for (int j = 0; j < Results.size(); j++) {
				//cout << "regionStart = " << Results[j].regionStart << "\tRegionEnd = " << Results[j].regionEnd << "\tName = " << Results[j].parent << "\tPerQP = " << Results[j].queryToParent << "\tLocalPerQP = " << Results[j].queryToParentLocal << "\tdivR = " << Results[j].divR << endl;
			//}
			
			//if (chimeraFlag == "yes") {
			
				//get sequence that were given from maligner results
				vector<SeqDist> seqs;
				map<string, float> removeDups;
				map<string, float>::iterator itDup;
				for (int j = 0; j < Results.size(); j++) {
					float dist = (Results[j].regionEnd - Results[j].regionStart + 1) * Results[j].queryToParentLocal;
					//only add if you are not a duplicate
					itDup = removeDups.find(Results[j].parent);
					if (itDup == removeDups.end()) { //this is not duplicate
						removeDups[Results[j].parent] = dist;
					}else if (dist > itDup->second) { //is this a stronger number for this parent
						removeDups[Results[j].parent] = dist;
					}
				}
				
				for (itDup = removeDups.begin(); itDup != removeDups.end(); itDup++) {
					Sequence* seq = getSequence(itDup->first); //makes copy so you can filter and mask and not effect template
					
					SeqDist member;
					member.seq = seq;
					member.dist = itDup->second;
					
					seqs.push_back(member);
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
//ofstream out;
//string name = querySeqs[i]->getName();
//cout << name << endl;
//name = name.substr(name.find_first_of("|")+1);
//cout << name << endl;
//name = name.substr(name.find_first_of("|")+1);
//cout << name << endl;
//name = name.substr(0, name.find_last_of("|"));
//cout << name << endl;
//string filename = name + ".seqsforslayer";
//openOutputFile(filename, out);	
//for (int u = 0; u < seqsForSlayer.size(); u++) { seqsForSlayer[u]->printSequence(out);	}
//out.close();
//filename = name + ".fasta";
//openOutputFile(filename, out);	
//q->printSequence(out);
//out.close();

			
				//mask then send to slayer...
				if (seqMask != "") {
					decalc->setMask(seqMask);

					//mask querys
					decalc->runMask(querySeqs[i]);
					
					//mask parents
					for (int k = 0; k < seqsForSlayer.size(); k++) {
						decalc->runMask(seqsForSlayer[k]);
					}
					
					for (int i = 0; i < numSeqs; i++) {
						spotMap[i] = decalc->getMaskMap();
					}

				}
				
				//send to slayer
				chimeraFlags[i] = slayer->getResults(querySeqs[i], seqsForSlayer);
				chimeraResults[i] = slayer->getOutput();
			
				//free memory
				for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }
			//}
			
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
void ChimeraSlayer::printBlock(data_struct data, ostream& out, int i){
	try {
		
		out << "parentA: " << data.parentA.getName() << "  parentB: " << data.parentB.getName()  << endl;
		out << "Left Window: " << spotMap[i][data.winLStart] << " " << spotMap[i][data.winLEnd] << endl;
		out << "Right Window: " << spotMap[i][data.winRStart] << " " << spotMap[i][data.winREnd] << endl;
		
		out << "Divergence of Query to Leftside ParentA and Rightside ParentB: " << data.divr_qla_qrb << '\t' << "Bootstrap: " << data.bsa << endl;
		out << "Divergence of Query to Rightside ParentA and Leftside ParentB: " << data.divr_qlb_qra << '\t' << "Bootstrap: " << data.bsb << endl;
		
		out << "Similarity of parents: " << data.ab << endl;
		out << "Similarity of query to parentA: " << data.qa << endl;
		out << "Similarity of query to parentB: " << data.qb << endl;
		
		//out << "Per_id(QL,A): " << data.qla << endl;
		//out << "Per_id(QL,B): " << data.qlb << endl;
		//out << "Per_id(QR,A): " << data.qra << endl;
		//out << "Per_id(QR,B): " << data.qrb << endl;

		
		out << "DeltaL: " << (data.qla - data.qlb) << endl;
		out << "DeltaR: " << (data.qra - data.qrb) << endl;

	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "printBlock");
		exit(1);
	}
}
//***************************************************************************************************************

