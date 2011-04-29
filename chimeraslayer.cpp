/*
 *  chimeraslayer.cpp
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Pschloss Lab. All rights reserved.
 *
 */

#include "chimeraslayer.h"
#include "chimerarealigner.h"
#include "kmerdb.hpp"
#include "blastdb.hpp"

//***************************************************************************************************************
ChimeraSlayer::ChimeraSlayer(string file, string temp, bool trim, string mode, int k, int ms, int mms, int win, float div, 
int minsim, int mincov, int minbs, int minsnp, int par, int it, int inc, int numw, bool r) : Chimera()  {  	
	try {
		fastafile = file;
		templateFileName = temp; templateSeqs = readSeqs(temp);
		searchMethod = mode;
		kmerSize = k;
		match = ms;
		misMatch = mms;
		window = win;
		divR = div;
		minSim = minsim;
		minCov = mincov;
		minBS = minbs;
		minSNP = minsnp;
		parents = par;
		iters = it;
		increment = inc;
		numWanted = numw;
		realign = r; 
		trimChimera = trim;
	
		decalc = new DeCalculator();	
		
		doPrep();
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "ChimeraSlayer");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraSlayer::ChimeraSlayer(string file, string temp, bool trim, map<string, int>& prior, string mode, int k, int ms, int mms, int win, float div, 
							 int minsim, int mincov, int minbs, int minsnp, int par, int it, int inc, int numw, bool r) : Chimera()  {  	
	try {
		fastafile = file; templateSeqs = readSeqs(fastafile);
		templateFileName = temp; 
		searchMethod = mode;
		kmerSize = k;
		match = ms;
		misMatch = mms;
		window = win;
		divR = div;
		minSim = minsim;
		minCov = mincov;
		minBS = minbs;
		minSNP = minsnp;
		parents = par;
		iters = it;
		increment = inc;
		numWanted = numw;
		realign = r; 
		trimChimera = trim;
		priority = prior;
		
		decalc = new DeCalculator();	
		
		createFilter(templateSeqs, 0.0); //just removed columns where all seqs have a gap
		
		//run filter on template
		for (int i = 0; i < templateSeqs.size(); i++) {  if (m->control_pressed) {  break; }  runFilter(templateSeqs[i]);  }

		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "ChimeraSlayer");
		exit(1);
	}
}
//***************************************************************************************************************
int ChimeraSlayer::doPrep() {
	try {
		if (searchMethod == "distance") { 
			//read in all query seqs
			vector<Sequence*> tempQuerySeqs = readSeqs(fastafile);
		
			vector<Sequence*> temp = templateSeqs;
			for (int i = 0; i < tempQuerySeqs.size(); i++) {  temp.push_back(tempQuerySeqs[i]);  }
		
			createFilter(temp, 0.0); //just removed columns where all seqs have a gap
		
			for (int i = 0; i < tempQuerySeqs.size(); i++) { delete tempQuerySeqs[i];  }
		
			if (m->control_pressed) {  return 0; } 
		
			//run filter on template copying templateSeqs into filteredTemplateSeqs
			for (int i = 0; i < templateSeqs.size(); i++) {  
				if (m->control_pressed) {  return 0; }
				
				Sequence* newSeq = new Sequence(templateSeqs[i]->getName(), templateSeqs[i]->getAligned());
				filteredTemplateSeqs.push_back(newSeq);
				runFilter(newSeq);  
			}
		}
		string 	kmerDBNameLeft;
		string 	kmerDBNameRight;
	
		//generate the kmerdb to pass to maligner
		if (searchMethod == "kmer") { 
			string templatePath = m->hasPath(templateFileName);
			string rightTemplateFileName = templatePath + "right." + m->getRootName(m->getSimpleName(templateFileName));
			databaseRight = new KmerDB(rightTemplateFileName, kmerSize);
				
			string leftTemplateFileName = templatePath + "left." + m->getRootName(m->getSimpleName(templateFileName));
			databaseLeft = new KmerDB(leftTemplateFileName, kmerSize);	
		#ifdef USE_MPI
			for (int i = 0; i < templateSeqs.size(); i++) {
					
				if (m->control_pressed) { return 0; } 
					
				string leftFrag = templateSeqs[i]->getUnaligned();
				leftFrag = leftFrag.substr(0, int(leftFrag.length() * 0.33));
					
				Sequence leftTemp(templateSeqs[i]->getName(), leftFrag);
				databaseLeft->addSequence(leftTemp);	
			}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(templateSeqs.size());
			
			for (int i = 0; i < templateSeqs.size(); i++) {
				if (m->control_pressed) { return 0; } 
					
				string rightFrag = templateSeqs[i]->getUnaligned();
				rightFrag = rightFrag.substr(int(rightFrag.length() * 0.66));
					
				Sequence rightTemp(templateSeqs[i]->getName(), rightFrag);
				databaseRight->addSequence(rightTemp);	
			}
			databaseRight->generateDB();
			databaseRight->setNumSeqs(templateSeqs.size());

		#else	
			//leftside
			kmerDBNameLeft = leftTemplateFileName.substr(0,leftTemplateFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTestLeft(kmerDBNameLeft.c_str());
			bool needToGenerateLeft = true;
			
			if(kmerFileTestLeft){	
				bool GoodFile = m->checkReleaseVersion(kmerFileTestLeft, m->getVersion());
				if (GoodFile) {  needToGenerateLeft = false;	}
			}
			
			if(needToGenerateLeft){	
			
				for (int i = 0; i < templateSeqs.size(); i++) {
					
					if (m->control_pressed) { return 0; } 
					
					string leftFrag = templateSeqs[i]->getUnaligned();
					leftFrag = leftFrag.substr(0, int(leftFrag.length() * 0.33));
					
					Sequence leftTemp(templateSeqs[i]->getName(), leftFrag);
					databaseLeft->addSequence(leftTemp);	
				}
				databaseLeft->generateDB();
				
			}else {	
				databaseLeft->readKmerDB(kmerFileTestLeft);	
			}
			kmerFileTestLeft.close();
			
			databaseLeft->setNumSeqs(templateSeqs.size());
			
			//rightside
			kmerDBNameRight = rightTemplateFileName.substr(0,rightTemplateFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTestRight(kmerDBNameRight.c_str());
			bool needToGenerateRight = true;
			
			if(kmerFileTestRight){	
				bool GoodFile = m->checkReleaseVersion(kmerFileTestRight, m->getVersion());
				if (GoodFile) {  needToGenerateRight = false;	}
			}
			
			if(needToGenerateRight){	
			
				for (int i = 0; i < templateSeqs.size(); i++) {
					if (m->control_pressed) { return 0; } 
					
					string rightFrag = templateSeqs[i]->getUnaligned();
					rightFrag = rightFrag.substr(int(rightFrag.length() * 0.66));
					
					Sequence rightTemp(templateSeqs[i]->getName(), rightFrag);
					databaseRight->addSequence(rightTemp);	
				}
				databaseRight->generateDB();
				
			}else {	
				databaseRight->readKmerDB(kmerFileTestRight);	
			}
			kmerFileTestRight.close();
			
			databaseRight->setNumSeqs(templateSeqs.size());
		#endif	
		}else if (searchMethod == "blast") {
		
			//generate blastdb
			databaseLeft = new BlastDB(-1.0, -1.0, 1, -3);

			for (int i = 0; i < templateSeqs.size(); i++) { 	databaseLeft->addSequence(*templateSeqs[i]);	}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(templateSeqs.size());
		}
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "doprep");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence*> ChimeraSlayer::getTemplate(Sequence* q, vector<Sequence*>& userTemplateFiltered) {
	try {
		
		//when template=self, the query file is sorted from most abundance to least abundant
		//userTemplate grows as the query file is processed by adding sequences that are not chimeric and more abundant
		vector<Sequence*> userTemplate;
		
		int myAbund = priority[q->getName()];
		
		for (int i = 0; i < templateSeqs.size(); i++) {
			
			if (m->control_pressed) { return userTemplate; } 
			
			//have I reached a sequence with the same abundance as myself?
			if (!(priority[templateSeqs[i]->getName()] > myAbund)) { break; }
			
			//if its am not chimeric add it
			if (chimericSeqs.count(templateSeqs[i]->getName()) == 0) { 
				userTemplate.push_back(templateSeqs[i]); 
				if (searchMethod == "distance") { userTemplateFiltered.push_back(filteredTemplateSeqs[i]); }
			}
		}
		
		string 	kmerDBNameLeft;
		string 	kmerDBNameRight;
		
		//generate the kmerdb to pass to maligner
		if (searchMethod == "kmer") { 
			string templatePath = m->hasPath(templateFileName);
			string rightTemplateFileName = templatePath + "right." + m->getRootName(m->getSimpleName(templateFileName));
			databaseRight = new KmerDB(rightTemplateFileName, kmerSize);
			
			string leftTemplateFileName = templatePath + "left." + m->getRootName(m->getSimpleName(templateFileName));
			databaseLeft = new KmerDB(leftTemplateFileName, kmerSize);	
#ifdef USE_MPI
			for (int i = 0; i < userTemplate.size(); i++) {
				
				if (m->control_pressed) { return userTemplate; } 
				
				string leftFrag = userTemplate[i]->getUnaligned();
				leftFrag = leftFrag.substr(0, int(leftFrag.length() * 0.33));
				
				Sequence leftTemp(userTemplate[i]->getName(), leftFrag);
				databaseLeft->addSequence(leftTemp);	
			}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(userTemplate.size());
			
			for (int i = 0; i < userTemplate.size(); i++) {
				if (m->control_pressed) { return userTemplate; } 
				
				string rightFrag = userTemplate[i]->getUnaligned();
				rightFrag = rightFrag.substr(int(rightFrag.length() * 0.66));
				
				Sequence rightTemp(userTemplate[i]->getName(), rightFrag);
				databaseRight->addSequence(rightTemp);	
			}
			databaseRight->generateDB();
			databaseRight->setNumSeqs(userTemplate.size());
			
#else	
			
			
			for (int i = 0; i < userTemplate.size(); i++) {
				
				if (m->control_pressed) { return userTemplate; } 
				
				string leftFrag = userTemplate[i]->getUnaligned();
				leftFrag = leftFrag.substr(0, int(leftFrag.length() * 0.33));
				
				Sequence leftTemp(userTemplate[i]->getName(), leftFrag);
				databaseLeft->addSequence(leftTemp);	
			}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(userTemplate.size());
				
			for (int i = 0; i < userTemplate.size(); i++) {
				if (m->control_pressed) { return userTemplate; }  
					
				string rightFrag = userTemplate[i]->getUnaligned();
				rightFrag = rightFrag.substr(int(rightFrag.length() * 0.66));
					
				Sequence rightTemp(userTemplate[i]->getName(), rightFrag);
				databaseRight->addSequence(rightTemp);	
			}
			databaseRight->generateDB();
			databaseRight->setNumSeqs(userTemplate.size());
#endif	
		}else if (searchMethod == "blast") {
			
			//generate blastdb
			databaseLeft = new BlastDB(-1.0, -1.0, 1, -3);

			for (int i = 0; i < userTemplate.size(); i++) { if (m->control_pressed) { return userTemplate; }   databaseLeft->addSequence(*userTemplate[i]);	}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(userTemplate.size());
		}
		
		return userTemplate;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getTemplate");
		exit(1);
	}
}

//***************************************************************************************************************
ChimeraSlayer::~ChimeraSlayer() { 	
	delete decalc;  
	if (templateFileName != "self") {
		if (searchMethod == "kmer") {  delete databaseRight;  delete databaseLeft;  }	
		else if (searchMethod == "blast") {  delete databaseLeft; }
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printHeader(ostream& out) {
	m->mothurOutEndLine();
	m->mothurOut("Only reporting sequence supported by " + toString(minBS) + "% of bootstrapped results.");
	m->mothurOutEndLine();
	
	out << "Name\tLeftParent\tRightParent\tDivQLAQRB\tPerIDQLAQRB\tBootStrapA\tDivQLBQRA\tPerIDQLBQRA\tBootStrapB\tFlag\tLeftWindow\tRightWindow\n";
}
//***************************************************************************************************************
Sequence* ChimeraSlayer::print(ostream& out, ostream& outAcc) {
	try {
		Sequence* trim = NULL;
		if (trimChimera) { trim = new Sequence(trimQuery.getName(), trimQuery.getAligned()); }
		
		if (chimeraFlags == "yes") {
			string chimeraFlag = "no";
			if(  (chimeraResults[0].bsa >= minBS && chimeraResults[0].divr_qla_qrb >= divR)
			   ||
			   (chimeraResults[0].bsb >= minBS && chimeraResults[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			
			
			if (chimeraFlag == "yes") {	
				if ((chimeraResults[0].bsa >= minBS) || (chimeraResults[0].bsb >= minBS)) {
					m->mothurOut(querySeq->getName() + "\tyes"); m->mothurOutEndLine();
					outAcc << querySeq->getName() << endl;
					
					if (templateFileName == "self") {  chimericSeqs.insert(querySeq->getName()); }
					
					if (trimChimera) {  
						int lengthLeft = chimeraResults[0].winLEnd - chimeraResults[0].winLStart;
						int lengthRight = chimeraResults[0].winREnd - chimeraResults[0].winRStart;
						
						string newAligned = trim->getAligned();

						if (lengthLeft > lengthRight) { //trim right
							for (int i = (chimeraResults[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
						}else { //trim left
							for (int i = 0; i < chimeraResults[0].winLEnd; i++) { newAligned[i] = '.'; }
						}
						trim->setAligned(newAligned);
					}
				}
			}
			
			printBlock(chimeraResults[0], chimeraFlag, out);
			out << endl;
		}else {  
			out << querySeq->getName() << "\tno" << endl; 
		}
		
		return trim;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence* ChimeraSlayer::print(ostream& out, ostream& outAcc, data_results leftPiece, data_results rightPiece) {
	try {
		Sequence* trim = NULL;
				
		if (trimChimera) { 
			string aligned = leftPiece.trimQuery.getAligned() + rightPiece.trimQuery.getAligned();
			trim = new Sequence(leftPiece.trimQuery.getName(), aligned); 
		}
		
		if ((leftPiece.flag == "yes") || (rightPiece.flag == "yes")) {
			
			string chimeraFlag = "no";
			if (leftPiece.flag == "yes") {
				
				if(  (leftPiece.results[0].bsa >= minBS && leftPiece.results[0].divr_qla_qrb >= divR)
					||
					(leftPiece.results[0].bsb >= minBS && leftPiece.results[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			}
			
			if (rightPiece.flag == "yes") {
				if ( (rightPiece.results[0].bsa >= minBS && rightPiece.results[0].divr_qla_qrb >= divR)
				 ||
				 (rightPiece.results[0].bsb >= minBS && rightPiece.results[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			}
			
			bool rightChimeric = false;
			bool leftChimeric = false;
			
			if (chimeraFlag == "yes") {	
				//which peice is chimeric or are both
				if (rightPiece.flag == "yes") { if ((rightPiece.results[0].bsa >= minBS) || (rightPiece.results[0].bsb >= minBS)) { rightChimeric = true; } }
				if (leftPiece.flag == "yes") { if ((leftPiece.results[0].bsa >= minBS) || (leftPiece.results[0].bsb >= minBS))	{ leftChimeric = true;	} }
				
				if (rightChimeric || leftChimeric) {
					m->mothurOut(querySeq->getName() + "\tyes"); m->mothurOutEndLine();
					outAcc << querySeq->getName() << endl;
					
					if (templateFileName == "self") {  chimericSeqs.insert(querySeq->getName()); }
					
					if (trimChimera) {  
						string newAligned = trim->getAligned();
												
						//right side is fine so keep that
						if ((leftChimeric) && (!rightChimeric)) {
							for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
						}else if ((!leftChimeric) && (rightChimeric)) { //leftside is fine so keep that
							for (int i = (rightPiece.results[0].winLStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
						}else { //both sides are chimeric, keep longest piece
							
							int lengthLeftLeft = leftPiece.results[0].winLEnd - leftPiece.results[0].winLStart;
							int lengthLeftRight = leftPiece.results[0].winREnd - leftPiece.results[0].winRStart;
							
							int longest = 1; // leftleft = 1, leftright = 2, rightleft = 3 rightright = 4
							int length = lengthLeftLeft;
							if (lengthLeftLeft < lengthLeftRight) { longest = 2;  length = lengthLeftRight; }
							
							int lengthRightLeft = rightPiece.results[0].winLEnd - rightPiece.results[0].winLStart;
							int lengthRightRight = rightPiece.results[0].winREnd - rightPiece.results[0].winRStart;
							
							if (lengthRightLeft > length) { longest = 3; length = lengthRightLeft;  }
							if (lengthRightRight > length) { longest = 4; }
							
							if (longest == 1) { //leftleft
								for (int i = (leftPiece.results[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else if (longest == 2) { //leftright
								//get rid of leftleft
								for (int i = (leftPiece.results[0].winLStart-1); i < (leftPiece.results[0].winLEnd-1); i++) { newAligned[i] = '.'; }
								//get rid of right
								for (int i = (rightPiece.results[0].winLStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else if (longest == 3) { //rightleft
								//get rid of left
								for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
								//get rid of rightright
								for (int i = (rightPiece.results[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else { //rightright
								//get rid of left
								for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
								//get rid of rightleft
								for (int i = (rightPiece.results[0].winLStart-1); i < (rightPiece.results[0].winLEnd-1); i++) { newAligned[i] = '.'; }
							}
						}
							
						trim->setAligned(newAligned);
					}
					
				}
			}
			
			printBlock(leftPiece, rightPiece, leftChimeric, rightChimeric, chimeraFlag, out);
			out << endl;
		}else {  
			out << querySeq->getName() << "\tno" << endl;  
		}
		
		return trim;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}

#ifdef USE_MPI
//***************************************************************************************************************
Sequence* ChimeraSlayer::print(MPI_File& out, MPI_File& outAcc, data_results leftPiece, data_results rightPiece) {
	try {
		MPI_Status status;
		bool results = false;
		string outAccString = "";
		string outputString = "";
		
		Sequence* trim = NULL;
		
		if (trimChimera) { 
			string aligned = leftPiece.trimQuery.getAligned() + rightPiece.trimQuery.getAligned();
			trim = new Sequence(leftPiece.trimQuery.getName(), aligned); 
		}
		
		
		if ((leftPiece.flag == "yes") || (rightPiece.flag == "yes")) {
			
			string chimeraFlag = "no";
			if (leftPiece.flag == "yes") {
				
				if(  (leftPiece.results[0].bsa >= minBS && leftPiece.results[0].divr_qla_qrb >= divR)
				   ||
				   (leftPiece.results[0].bsb >= minBS && leftPiece.results[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			}
			
			if (rightPiece.flag == "yes") {
				if ( (rightPiece.results[0].bsa >= minBS && rightPiece.results[0].divr_qla_qrb >= divR)
					||
					(rightPiece.results[0].bsb >= minBS && rightPiece.results[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			}
			
			bool rightChimeric = false;
			bool leftChimeric = false;
			
			if (chimeraFlag == "yes") {	
				//which peice is chimeric or are both
				if (rightPiece.flag == "yes") { if ((rightPiece.results[0].bsa >= minBS) || (rightPiece.results[0].bsb >= minBS)) { rightChimeric = true; } }
				if (leftPiece.flag == "yes") { if ((leftPiece.results[0].bsa >= minBS) || (leftPiece.results[0].bsb >= minBS))	{ leftChimeric = true;	} }
				
				if (rightChimeric || leftChimeric) {
					cout << querySeq->getName() <<  "\tyes" << endl;
					outAccString += querySeq->getName() + "\n";
					results = true;
					
					if (templateFileName == "self") {  chimericSeqs.insert(querySeq->getName()); }
					
					//write to accnos file
					int length = outAccString.length();
					char* buf2 = new char[length];
					memcpy(buf2, outAccString.c_str(), length);
				
					MPI_File_write_shared(outAcc, buf2, length, MPI_CHAR, &status);
					delete buf2;
					
					if (trimChimera) {  
						string newAligned = trim->getAligned();
						
						//right side is fine so keep that
						if ((leftChimeric) && (!rightChimeric)) {
							for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
						}else if ((!leftChimeric) && (rightChimeric)) { //leftside is fine so keep that
							for (int i = (rightPiece.results[0].winLStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
						}else { //both sides are chimeric, keep longest piece
							
							int lengthLeftLeft = leftPiece.results[0].winLEnd - leftPiece.results[0].winLStart;
							int lengthLeftRight = leftPiece.results[0].winREnd - leftPiece.results[0].winRStart;
							
							int longest = 1; // leftleft = 1, leftright = 2, rightleft = 3 rightright = 4
							int length = lengthLeftLeft;
							if (lengthLeftLeft < lengthLeftRight) { longest = 2;  length = lengthLeftRight; }
							
							int lengthRightLeft = rightPiece.results[0].winLEnd - rightPiece.results[0].winLStart;
							int lengthRightRight = rightPiece.results[0].winREnd - rightPiece.results[0].winRStart;
							
							if (lengthRightLeft > length) { longest = 3; length = lengthRightLeft;  }
							if (lengthRightRight > length) { longest = 4; }
							
							if (longest == 1) { //leftleft
								for (int i = (leftPiece.results[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else if (longest == 2) { //leftright
								//get rid of leftleft
								for (int i = (leftPiece.results[0].winLStart-1); i < (leftPiece.results[0].winLEnd-1); i++) { newAligned[i] = '.'; }
								//get rid of right
								for (int i = (rightPiece.results[0].winLStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else if (longest == 3) { //rightleft
								//get rid of left
								for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
								//get rid of rightright
								for (int i = (rightPiece.results[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else { //rightright
								//get rid of left
								for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
								//get rid of rightleft
								for (int i = (rightPiece.results[0].winLStart-1); i < (rightPiece.results[0].winLEnd-1); i++) { newAligned[i] = '.'; }
							}
						}
						
						trim->setAligned(newAligned);
					}
					
				}
			}
			
			outputString = getBlock(leftPiece, rightPiece, leftChimeric, rightChimeric, chimeraFlag);
			outputString += "\n";
		
			//write to output file
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
				
			MPI_File_write_shared(out, buf, length, MPI_CHAR, &status);
			delete buf;

		}else {  
			outputString += querySeq->getName() + "\tno\n";  
	
			//write to output file
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
				
			MPI_File_write_shared(out, buf, length, MPI_CHAR, &status);
			delete buf;
		}
		
		
		return trim;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence* ChimeraSlayer::print(MPI_File& out, MPI_File& outAcc) {
	try {
		MPI_Status status;
		bool results = false;
		string outAccString = "";
		string outputString = "";
		
		Sequence* trim = NULL;
		if (trimChimera) { trim = new Sequence(trimQuery.getName(), trimQuery.getAligned()); }
		
		if (chimeraFlags == "yes") {
			string chimeraFlag = "no";
			if(  (chimeraResults[0].bsa >= minBS && chimeraResults[0].divr_qla_qrb >= divR)
			   ||
			   (chimeraResults[0].bsb >= minBS && chimeraResults[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			
			
			if (chimeraFlag == "yes") {	
				if ((chimeraResults[0].bsa >= minBS) || (chimeraResults[0].bsb >= minBS)) {
					cout << querySeq->getName() <<  "\tyes" << endl;
					outAccString += querySeq->getName() + "\n";
					results = true;
					
					if (templateFileName == "self") {  chimericSeqs.insert(querySeq->getName()); }
					
					//write to accnos file
					int length = outAccString.length();
					char* buf2 = new char[length];
					memcpy(buf2, outAccString.c_str(), length);
					
					MPI_File_write_shared(outAcc, buf2, length, MPI_CHAR, &status);
					delete buf2;
					
					if (trimChimera) {  
						int lengthLeft = chimeraResults[0].winLEnd - chimeraResults[0].winLStart;
						int lengthRight = chimeraResults[0].winREnd - chimeraResults[0].winRStart;
						
						string newAligned = trim->getAligned();
						if (lengthLeft > lengthRight) { //trim right
							for (int i = (chimeraResults[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
						}else { //trim left
							for (int i = 0; i < (chimeraResults[0].winLEnd-1); i++) { newAligned[i] = '.'; }
						}
						trim->setAligned(newAligned);	
					}
				}
			}
			
			outputString = getBlock(chimeraResults[0], chimeraFlag);
			outputString += "\n";
			
			//write to output file
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write_shared(out, buf, length, MPI_CHAR, &status);
			delete buf;
			
		}else {  
			outputString += querySeq->getName() + "\tno\n";  
			
			//write to output file
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write_shared(out, buf, length, MPI_CHAR, &status);
			delete buf;
		}
		
		return trim;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}
#endif

//***************************************************************************************************************
int ChimeraSlayer::getChimeras(Sequence* query) {
	try {
		
		trimQuery.setName(query->getName()); trimQuery.setAligned(query->getAligned());
		printResults.trimQuery = trimQuery; 
		
		chimeraFlags = "no";
		printResults.flag = "no";
		
		querySeq = query;
		
		//you must create a template
		vector<Sequence*> thisTemplate;
		vector<Sequence*> thisFilteredTemplate;
		if (templateFileName != "self") { thisTemplate = templateSeqs; thisFilteredTemplate = filteredTemplateSeqs; }
		else {  thisTemplate = getTemplate(query, thisFilteredTemplate);  } //fills this template and creates the databases
		
		if (m->control_pressed) {  return 0;  }
		
		if (thisTemplate.size() == 0) {  return 0; } //not chimeric
		
		//moved this out of maligner - 4/29/11
		vector<Sequence*> refSeqs = getRefSeqs(query, thisTemplate, thisFilteredTemplate);
			
		Maligner maligner(refSeqs, match, misMatch, divR, minSim, minCov); 
		Slayer slayer(window, increment, minSim, divR, iters, minSNP, minBS);
		
		if (templateFileName == "self") {
			if (searchMethod == "kmer") {  delete databaseRight;  delete databaseLeft;  }	
			else if (searchMethod == "blast") {  delete databaseLeft; }
		}
	
		if (m->control_pressed) {  return 0;  }

		string chimeraFlag = maligner.getResults(query, decalc);
		
		if (m->control_pressed) {  return 0;  }
		
		vector<results> Results = maligner.getOutput();
		
		for (int i = 0; i < refSeqs.size(); i++) {  delete refSeqs[i];	}
		
		if (chimeraFlag == "yes") {
			
			if (realign) {
				ChimeraReAligner realigner(thisTemplate, match, misMatch);
				realigner.reAlign(query, Results);
			}
			
			//get sequence that were given from maligner results
			vector<SeqDist> seqs;
			map<string, float> removeDups;
			map<string, float>::iterator itDup;
			map<string, string> parentNameSeq;
			map<string, string>::iterator itSeq;
			for (int j = 0; j < Results.size(); j++) {
				float dist = (Results[j].regionEnd - Results[j].regionStart + 1) * Results[j].queryToParentLocal;
				//only add if you are not a duplicate
				itDup = removeDups.find(Results[j].parent);
				if (itDup == removeDups.end()) { //this is not duplicate
					removeDups[Results[j].parent] = dist;
					parentNameSeq[Results[j].parent] = Results[j].parentAligned;
				}else if (dist > itDup->second) { //is this a stronger number for this parent
					removeDups[Results[j].parent] = dist;
					parentNameSeq[Results[j].parent] = Results[j].parentAligned;
				}
			}
			
			for (itDup = removeDups.begin(); itDup != removeDups.end(); itDup++) {
				itSeq = parentNameSeq.find(itDup->first);
				Sequence* seq = new Sequence(itDup->first, itSeq->second);
				
				SeqDist member;
				member.seq = seq;
				member.dist = itDup->second;
				
				seqs.push_back(member);
			}
			
			//limit number of parents to explore - default 3
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
			
			if (m->control_pressed) {  for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }  return 0;  }

			//send to slayer
			chimeraFlags = slayer.getResults(query, seqsForSlayer);
			if (m->control_pressed) {  return 0;  }
			chimeraResults = slayer.getOutput();
			
			printResults.flag = chimeraFlags;
			printResults.results = chimeraResults;
			
			//free memory
			for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printBlock(data_struct data, string flag, ostream& out){
	try {
		out << querySeq->getName() << '\t';
		out << data.parentA.getName() << "\t" << data.parentB.getName()  << '\t';
	
		out << data.divr_qla_qrb << '\t' << data.qla_qrb << '\t' << data.bsa << '\t';
		out << data.divr_qlb_qra << '\t' << data.qlb_qra << '\t' << data.bsb << '\t';
		
		out << flag << '\t' << data.winLStart << "-" << data.winLEnd << '\t' << data.winRStart << "-" << data.winREnd << '\t';
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "printBlock");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printBlock(data_results leftdata, data_results rightdata, bool leftChimeric, bool rightChimeric, string flag, ostream& out){
	try {
		
		if ((leftChimeric) && (!rightChimeric)) { //print left
			out << querySeq->getName() << '\t';
			out << leftdata.results[0].parentA.getName() << "\t" << leftdata.results[0].parentB.getName()  << '\t';
			
			out << leftdata.results[0].divr_qla_qrb << '\t' << leftdata.results[0].qla_qrb << '\t' << leftdata.results[0].bsa << '\t';
			out << leftdata.results[0].divr_qlb_qra << '\t' << leftdata.results[0].qlb_qra << '\t' << leftdata.results[0].bsb << '\t';
		
			out << flag << '\t' << leftdata.results[0].winLStart << "-" << leftdata.results[0].winLEnd << '\t' << leftdata.results[0].winRStart << "-" << leftdata.results[0].winREnd << '\t';
		
		}else if ((!leftChimeric) && (rightChimeric)) {  //print right
			out << querySeq->getName() << '\t';
			out << rightdata.results[0].parentA.getName() << "\t" << rightdata.results[0].parentB.getName()  << '\t';
			
			out << rightdata.results[0].divr_qla_qrb << '\t' << rightdata.results[0].qla_qrb << '\t' << rightdata.results[0].bsa << '\t';
			out << rightdata.results[0].divr_qlb_qra << '\t' << rightdata.results[0].qlb_qra << '\t' << rightdata.results[0].bsb << '\t';
			
			out << flag << '\t' << rightdata.results[0].winLStart << "-" << rightdata.results[0].winLEnd << '\t' << rightdata.results[0].winRStart << "-" << rightdata.results[0].winREnd << '\t';			
			
		}else  { //print both results
			if (leftdata.flag == "yes") {
				out << querySeq->getName() + "_LEFT" << '\t';
				out << leftdata.results[0].parentA.getName() << "\t" << leftdata.results[0].parentB.getName()  << '\t';
				
				out << leftdata.results[0].divr_qla_qrb << '\t' << leftdata.results[0].qla_qrb << '\t' << leftdata.results[0].bsa << '\t';
				out << leftdata.results[0].divr_qlb_qra << '\t' << leftdata.results[0].qlb_qra << '\t' << leftdata.results[0].bsb << '\t';
				
				out << flag << '\t' << leftdata.results[0].winLStart << "-" << leftdata.results[0].winLEnd << '\t' << leftdata.results[0].winRStart << "-" << leftdata.results[0].winREnd << '\t';
			}
			
			if (rightdata.flag == "yes") {
				if (leftdata.flag == "yes") { out << endl; }
				
				out << querySeq->getName() + "_RIGHT"<< '\t';
				out << rightdata.results[0].parentA.getName() << "\t" << rightdata.results[0].parentB.getName()  << '\t';
				
				out << rightdata.results[0].divr_qla_qrb << '\t' << rightdata.results[0].qla_qrb << '\t' << rightdata.results[0].bsa << '\t';
				out << rightdata.results[0].divr_qlb_qra << '\t' << rightdata.results[0].qlb_qra << '\t' << rightdata.results[0].bsb << '\t';
				
				out << flag << '\t' << rightdata.results[0].winLStart << "-" << rightdata.results[0].winLEnd << '\t' << rightdata.results[0].winRStart << "-" << rightdata.results[0].winREnd << '\t';			
		
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "printBlock");
		exit(1);
	}
}
//***************************************************************************************************************
string ChimeraSlayer::getBlock(data_results leftdata, data_results rightdata, bool leftChimeric, bool rightChimeric, string flag){
	try {
		
		string out = "";
		
		if ((leftChimeric) && (!rightChimeric)) { //get left
			out += querySeq->getName() + "\t";
			out += leftdata.results[0].parentA.getName() + "\t" + leftdata.results[0].parentB.getName() + "\t";
			
			out += toString(leftdata.results[0].divr_qla_qrb) + "\t" + toString(leftdata.results[0].qla_qrb) + "\t" + toString(leftdata.results[0].bsa) + "\t";
			out += toString(leftdata.results[0].divr_qlb_qra) + "\t" + toString(leftdata.results[0].qlb_qra) + "\t" + toString(leftdata.results[0].bsb) + "\t";
			
			out += flag + "\t" + toString(leftdata.results[0].winLStart) + "-" + toString(leftdata.results[0].winLEnd) + "\t" + toString(leftdata.results[0].winRStart) + "-" + toString(leftdata.results[0].winREnd) + "\t";
			
		}else if ((!leftChimeric) && (rightChimeric)) {  //print right
			out += querySeq->getName() + "\t";
			out += rightdata.results[0].parentA.getName() + "\t" + rightdata.results[0].parentB.getName()  + "\t";
			
			out += toString(rightdata.results[0].divr_qla_qrb) + "\t" + toString(rightdata.results[0].qla_qrb) + "\t" + toString(rightdata.results[0].bsa) + "\t";
			out += toString(rightdata.results[0].divr_qlb_qra) + "\t" + toString(rightdata.results[0].qlb_qra) + "\t" + toString(rightdata.results[0].bsb) + "\t";
			
			out += flag + "\t" + toString(rightdata.results[0].winLStart) + "-" + toString(rightdata.results[0].winLEnd) + "\t" + toString(rightdata.results[0].winRStart) + "-" + toString(rightdata.results[0].winREnd) + "\t";			
			
		}else  { //print both results
			
			if (leftdata.flag == "yes") {
				out += querySeq->getName() + "_LEFT\t";
				out += leftdata.results[0].parentA.getName() + "\t" + leftdata.results[0].parentB.getName() + "\t";
				
				out += toString(leftdata.results[0].divr_qla_qrb) + "\t" + toString(leftdata.results[0].qla_qrb) + "\t" + toString(leftdata.results[0].bsa) + "\t";
				out += toString(leftdata.results[0].divr_qlb_qra) + "\t" + toString(leftdata.results[0].qlb_qra) + "\t" + toString(leftdata.results[0].bsb) + "\t";
				
				out += flag + "\t" + toString(leftdata.results[0].winLStart) + "-" + toString(leftdata.results[0].winLEnd) + "\t" + toString(leftdata.results[0].winRStart) + "-" + toString(leftdata.results[0].winREnd) + "\t";
			}
			
			if (rightdata.flag == "yes") {
				if (leftdata.flag == "yes") { out += "\n"; }
				out +=  querySeq->getName() + "_RIGHT\t";
				out += rightdata.results[0].parentA.getName() + "\t" + rightdata.results[0].parentB.getName()  + "\t";
				
				out += toString(rightdata.results[0].divr_qla_qrb) + "\t" + toString(rightdata.results[0].qla_qrb) + "\t" + toString(rightdata.results[0].bsa) + "\t";
				out += toString(rightdata.results[0].divr_qlb_qra) + "\t" + toString(rightdata.results[0].qlb_qra) + "\t" + toString(rightdata.results[0].bsb) + "\t";
				
				out += flag + "\t" + toString(rightdata.results[0].winLStart) + "-" + toString(rightdata.results[0].winLEnd) + "\t" + toString(rightdata.results[0].winRStart) + "-" + toString(rightdata.results[0].winREnd) + "\t";			
			}
		}
		
		return out;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getBlock");
		exit(1);
	}
}
//***************************************************************************************************************
string ChimeraSlayer::getBlock(data_struct data, string flag){
	try {
		
		string outputString = "";
		
		outputString += querySeq->getName() + "\t";
		outputString += data.parentA.getName() + "\t" + data.parentB.getName()  + "\t";
			
		outputString += toString(data.divr_qla_qrb) + "\t" + toString(data.qla_qrb) + "\t" + toString(data.bsa) + "\t";
		outputString += toString(data.divr_qlb_qra) + "\t" + toString(data.qlb_qra) + "\t" + toString(data.bsb) + "\t";
		
		outputString += flag + "\t" + toString(data.winLStart) + "-" + toString(data.winLEnd) + "\t" + toString(data.winRStart) + "-" + toString(data.winREnd) + "\t";
		
		return outputString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getBlock");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence*> ChimeraSlayer::getRefSeqs(Sequence* q, vector<Sequence*>& thisTemplate, vector<Sequence*>& thisFilteredTemplate){
	try {
		
		vector<Sequence*> refSeqs;
		
		if (searchMethod == "distance") {
			//find closest seqs to query in template - returns copies of seqs so trim does not destroy - remember to deallocate
			Sequence* newSeq = new Sequence(q->getName(), q->getAligned());
			runFilter(newSeq);
			refSeqs = decalc->findClosest(newSeq, thisTemplate, thisFilteredTemplate, numWanted);
		}else if (searchMethod == "blast")  {
			refSeqs = getBlastSeqs(q, thisTemplate, numWanted); //fills indexes
		}else if (searchMethod == "kmer") {
			refSeqs = getKmerSeqs(q, thisTemplate, numWanted); //fills indexes
		}else { m->mothurOut("not valid search."); exit(1);  } //should never get here
		
		return refSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getRefSeqs");
		exit(1);
	}
}
//***************************************************************************************************************/
vector<Sequence*> ChimeraSlayer::getBlastSeqs(Sequence* q, vector<Sequence*>& db, int num) {
	try {	
		
		vector<Sequence*> refResults;
		
		//get parts of query
		string queryUnAligned = q->getUnaligned();
		string leftQuery = queryUnAligned.substr(0, int(queryUnAligned.length() * 0.33)); //first 1/3 of the sequence
		string rightQuery = queryUnAligned.substr(int(queryUnAligned.length() * 0.66)); //last 1/3 of the sequence
		
		Sequence* queryLeft = new Sequence(q->getName()+"left", leftQuery);
		Sequence* queryRight = new Sequence(q->getName()+"right", rightQuery);
		
		vector<int> tempIndexesLeft = databaseLeft->findClosestMegaBlast(queryLeft, num+1);
		vector<int> tempIndexesRight = databaseLeft->findClosestMegaBlast(queryRight, num+1);
		
		vector<int> smaller;
		vector<int> larger;
		
		if (tempIndexesRight.size() < tempIndexesLeft.size()) { smaller = tempIndexesRight;  larger = tempIndexesLeft;  }
		else { smaller = tempIndexesLeft;  larger = tempIndexesRight;  } 
		
		//merge results		
		map<int, int> seen;
		map<int, int>::iterator it;
		vector<int> mergedResults;
		for (int i = 0; i < smaller.size(); i++) {
			//add left if you havent already
			it = seen.find(smaller[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(smaller[i]);
				seen[smaller[i]] = smaller[i];
			}
			
			//add right if you havent already
			it = seen.find(larger[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(larger[i]);
				seen[larger[i]] = larger[i];
			}
		}
		
		for (int i = smaller.size(); i < larger.size(); i++) {
			//add right if you havent already
			it = seen.find(larger[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(larger[i]);
				seen[larger[i]] = larger[i];
			}
		}
		//numWanted = mergedResults.size();

		//cout << q->getName() << " merged results size = " << mergedResults.size() << '\t' << "numwanted = " << numWanted <<  endl;		
		for (int i = 0; i < mergedResults.size(); i++) {
			//cout << db[mergedResults[i]]->getName()  << '\t' << mergedResults[i] << endl;	
			
			if (db[mergedResults[i]]->getName() != q->getName()) { 
				Sequence* temp = new Sequence(db[mergedResults[i]]->getName(), db[mergedResults[i]]->getAligned());
				refResults.push_back(temp);
				//cout << db[mergedResults[i]]->getName() << endl;
			}
			
			//cout << mergedResults[i] << endl;
		}
		//cout << "done " << q->getName()  << endl;		
		delete queryRight;
		delete queryLeft;
		
		return refResults;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getBlastSeqs");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence*> ChimeraSlayer::getKmerSeqs(Sequence* q, vector<Sequence*>& db, int num) {
	try {	
		
		//get parts of query
		string queryUnAligned = q->getUnaligned();
		string leftQuery = queryUnAligned.substr(0, int(queryUnAligned.length() * 0.33)); //first 1/3 of the sequence
		string rightQuery = queryUnAligned.substr(int(queryUnAligned.length() * 0.66)); //last 1/3 of the sequence
		
		Sequence* queryLeft = new Sequence(q->getName(), leftQuery);
		Sequence* queryRight = new Sequence(q->getName(), rightQuery);
		
		vector<int> tempIndexesLeft = databaseLeft->findClosestSequences(queryLeft, numWanted);
		vector<int> tempIndexesRight = databaseRight->findClosestSequences(queryRight, numWanted);
		
		//merge results		
		map<int, int> seen;
		map<int, int>::iterator it;
			vector<int> mergedResults;
		for (int i = 0; i < tempIndexesLeft.size(); i++) {
			//add left if you havent already
			it = seen.find(tempIndexesLeft[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesLeft[i]);
				seen[tempIndexesLeft[i]] = tempIndexesLeft[i];
			}
			
			//add right if you havent already
			it = seen.find(tempIndexesRight[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesRight[i]);
				seen[tempIndexesRight[i]] = tempIndexesRight[i];
			}
		}
		
		//numWanted = mergedResults.size();
			
		//cout << q->getName() << endl;		
		vector<Sequence*> refResults;
		for (int i = 0; i < mergedResults.size(); i++) {
			//cout << db[mergedResults[i]]->getName() << endl;	
			if (db[mergedResults[i]]->getName() != q->getName()) { 
				Sequence* temp = new Sequence(db[mergedResults[i]]->getName(), db[mergedResults[i]]->getAligned());
				refResults.push_back(temp);
			}
		}
		//cout << endl;		
		delete queryRight;
		delete queryLeft;
		
		return refResults;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getKmerSeqs");
		exit(1);
	}
}
//***************************************************************************************************************


