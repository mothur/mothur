/*
 *  ccode.cpp
 *  Mothur
 *
 *  Created by westcott on 8/24/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "ccode.h"
#include "ignoregaps.h"
#include "eachgapdist.h"


//***************************************************************************************************************
Ccode::Ccode(string filename, string temp) {  fastafile = filename;  templateFile = temp;  }
//***************************************************************************************************************

Ccode::~Ccode() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)		{  delete querySeqs[i];		}
		for (int i = 0; i < templateSeqs.size(); i++)	{  delete templateSeqs[i];	}
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "~Ccode");
		exit(1);
	}
}	
//***************************************************************************************************************
void Ccode::print(ostream& out) {
	try {
		
		mothurOutEndLine();
		
		string mapInfo = getRootName(fastafile) + "mapinfo";
		ofstream out2;
		openOutputFile(mapInfo, out2);
		
		out2 << "Place in masked, filtered and trimmed sequence\tPlace in original alignment" << endl;
		
		for (int j = 0; j < querySeqs.size(); j++) {
			out2 << querySeqs[j]->getName() << endl;
			for (it = spotMap[j].begin(); it!= spotMap[j].end(); it++) {
				out2 << it->first << '\t' << it->second << endl;
			}
		}
		out2.close();
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		out << "For full window mapping info refer to " << mapInfo << endl << endl;
		
		for (int i = 0; i < querySeqs.size(); i++) {
		
			out << querySeqs[i]->getName() << endl << endl << "Reference sequences used and distance to query:" << endl;
			
			for (int j = 0; j < closest[i].size(); j++) {
				out << setprecision(3) << closest[i][j].seq->getName() << '\t' << closest[i][j].dist << endl;
			}
			out << endl << endl;
		
			//for each window
			//window mapping info.
			out << "Mapping information: ";
			//you mask and did not filter
			if ((seqMask != "") && (!filter)) { out << "mask and trim."; }
				
			//you filtered and did not mask
			if ((seqMask == "") && (filter)) { out << "filter and trim."; }
				
			//you masked and filtered
			if ((seqMask != "") && (filter)) { out << "mask, filter and trim."; }
			
			out << endl << "Window\tStartPos\tEndPos" << endl;
			it = trim[i].begin();
			
			for (int k = 0; k < windows[i].size()-1; k++) {
				out << k+1 << '\t' << spotMap[i][windows[i][k]-it->first] << '\t' << spotMap[i][windows[i][k]-it->first+windowSizes[i]] << endl;
			}
			
			out << windows[i].size() << '\t' << spotMap[i][windows[i][windows[i].size()-1]-it->first] << '\t' << spotMap[i][it->second-it->first-1] << endl;
			out << endl;
			
			out << "Window\tAvgQ\t(sdQ)\tAvgR\t(sdR)\tRatio\tAnova" << endl;
			for (int k = 0; k < windows[i].size(); k++) {
				float ds = averageQuery[i][k] / averageRef[i][k]; 
				out << k+1 << '\t' << averageQuery[i][k] << '\t' << sdQuery[i][k] << '\t' << averageRef[i][k] << '\t'<< sdRef[i][k] << '\t' << ds << '\t' << anova[i][k] << endl;
			}
			out << endl;
			
			//varRef
			//varQuery
			/* F test for differences among variances.
			* varQuery is expected to be higher or similar than varRef */
			//float fs = varQuery[query][i] / varRef[query][i];	/* F-Snedecor, test for differences of variances */
			
			bool results = false;	
						
			//confidence limit, t - Student, anova
			out << "Window\tConfidenceLimit\tt-Student\tAnova" << endl;
			
			for (int k = 0; k < windows[i].size(); k++) {
				string temp = "";
				if (isChimericConfidence[i][k]) {  temp += "*\t"; }
				else { temp += "\t"; }
				
				if (isChimericTStudent[i][k]) {  temp += "*\t"; }
				else { temp += "\t"; }
				
				if (isChimericANOVA[i][k]) {  temp += "*\t"; }
				else { temp += "\t"; }
			
				out << k+1 << '\t' << temp << endl;
				
				if (temp == "*\t*\t*\t") {  results = true;  }
			}
			out << endl;	
			
			if (results) {
				mothurOut(querySeqs[i]->getName() + " was found have at least one chimeric window."); mothurOutEndLine();
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "print");
		exit(1);
	}
}

//***************************************************************************************************************
void Ccode::getChimeras() {
	try {
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		templateSeqs = readSeqs(templateFile);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		closest.resize(numSeqs);
		
		refCombo.resize(numSeqs, 0);
		sumRef.resize(numSeqs); 
		varRef.resize(numSeqs); 
		varQuery.resize(numSeqs); 
		sdRef.resize(numSeqs); 
		sdQuery.resize(numSeqs);     
		sumQuery.resize(numSeqs);
		sumSquaredRef.resize(numSeqs); 
		sumSquaredQuery.resize(numSeqs); 
		averageRef.resize(numSeqs);
		averageQuery.resize(numSeqs);
		anova.resize(numSeqs);
		isChimericConfidence.resize(numSeqs);
		isChimericTStudent.resize(numSeqs);
		isChimericANOVA.resize(numSeqs);
		trim.resize(numSeqs);
		spotMap.resize(numSeqs);
		windowSizes.resize(numSeqs, window);
		windows.resize(numSeqs);
		
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
	
		distCalc = new eachGapDist();
		decalc = new DeCalculator();
		
		//find closest
		if (processors == 1) { 
			mothurOut("Finding top matches for sequences... "); cout.flush();
			closest = findClosest(lines[0]->start, lines[0]->end, numWanted);
			mothurOut("Done."); mothurOutEndLine();
		}else {		createProcessesClosest();		}

		//initialize spotMap
		for (int j = 0; j < numSeqs; j++) {
			for (int i = 0; i < querySeqs[0]->getAligned().length(); i++) {
				spotMap[j][i] = i;
			}
		}
		
		//mask sequences if the user wants to 
		if (seqMask != "") {
			decalc->setMask(seqMask);
			
			//mask querys
			for (int i = 0; i < querySeqs.size(); i++) {
				decalc->runMask(querySeqs[i]);
			}
		
			//mask templates
			for (int i = 0; i < templateSeqs.size(); i++) {
				decalc->runMask(templateSeqs[i]);
			}
			
			for (int i = 0; i < numSeqs; i++) {
				spotMap[i] = decalc->getMaskMap();
			}
		}
		
		if (filter) {
			vector<Sequence*> temp = templateSeqs;
			for (int i = 0; i < querySeqs.size(); i++) { temp.push_back(querySeqs[i]);  }
			
			createFilter(temp);
		
			runFilter(querySeqs);
			runFilter(templateSeqs);
			
			//update spotMap
			map<int, int> newMap;
			int spot = 0;
			int j = 0;
			
			for (int i = 0; i < filterString.length(); i++) {
				if (filterString[i] == '1') {
					//add to newMap
					newMap[spot] = spotMap[j][i];
					spot++;  
				}
			}
			
			for (int i = 0; i < numSeqs; i++) {
				spotMap[i] = newMap;
			}
		}

		//trim sequences - this follows ccodes remove_extra_gaps 
		for (int i = 0; i < querySeqs.size(); i++) {
			trimSequences(i);
		}
		
		//windows are equivalent to words - ccode paper recommends windows are between 5% and 20% on alignment length().  
		//Our default will be 10% and we will warn if user tries to use a window above or below these recommendations
		for (int i = 0; i < querySeqs.size(); i++) {
			windows[i] = findWindows(i);  
		}

		//remove sequences that are more than 20% different and less than 0.5% different - may want to allow user to specify this later 
		if (processors == 1) { 
			for (int i = 0; i < closest.size(); i++) {
				removeBadReferenceSeqs(closest[i], i);
			}
		}else {		createProcessesRemoveBad();		}

		
		if (processors == 1) { 
			//find the averages for each querys references
			for (int i = 0; i < numSeqs; i++) {
				getAverageRef(closest[i], i);  //fills sumRef[i], averageRef[i], sumSquaredRef[i] and refCombo[i].
			}
			
			//find the averages for the query 
			for (int i = 0; i < numSeqs; i++) {
				getAverageQuery(closest[i], i);  //fills sumQuery[i], averageQuery[i], sumSquaredQuery[i].
			}
		}else {		createProcessesAverages();		}
		
		if (processors == 1) { 
			//find the averages for each querys references 
			for (int i = 0; i < numSeqs; i++) {
				findVarianceRef(i);  //fills varRef[i] and sdRef[i] also sets minimum error rate to 0.001 to avoid divide by 0.
			}
			
			//find the averages for the query 
			for (int i = 0; i < numSeqs; i++) {
				findVarianceQuery(i);  //fills varQuery[i] and sdQuery[i] also sets minimum error rate to 0.001 to avoid divide by 0.
			}
		}else {		createProcessesVariances();		}
		
		if (processors == 1) { 
			for (int i = 0; i < numSeqs; i++) {
				determineChimeras(i);  //fills anova, isChimericConfidence[i], isChimericTStudent[i] and isChimericANOVA[i]. 
			}
		}else {		createProcessesDetermine();		}
		
		//free memory
		for (int i = 0; i < lines.size(); i++)					{	delete lines[i];				}
		delete distCalc;
		delete decalc;
			
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getChimeras");
		exit(1);
	}
}
/***************************************************************************************************************/
//ccode algo says it does this to "Removes the initial and final gaps to avoid biases due to incomplete sequences."
void Ccode::trimSequences(int query) {
	try {
		
		int frontPos = 0;  //should contain first position in all seqs that is not a gap character
		int rearPos = querySeqs[query]->getAligned().length();
		
		//********find first position in closest seqs that is a non gap character***********//
		//find first position all query seqs that is a non gap character
		for (int i = 0; i < closest[query].size(); i++) {
			
			string aligned = closest[query][i].seq->getAligned();
			int pos = 0;
			
			//find first spot in this seq
			for (int j = 0; j < aligned.length(); j++) {
				if (isalpha(aligned[j])) {
					pos = j;
					break;
				}
			}
			
			//save this spot if it is the farthest
			if (pos > frontPos) { frontPos = pos; }
		}
		
		//find first position all querySeq[query] that is a non gap character
		string aligned = querySeqs[query]->getAligned();
		int pos = 0;
			
		//find first spot in this seq
		for (int j = 0; j < aligned.length(); j++) {
			if (isalpha(aligned[j])) {
				pos = j;
				break;
			}
		}
		
		//save this spot if it is the farthest
		if (pos > frontPos) { frontPos = pos; }
		
		
		//********find last position in closest seqs that is a non gap character***********//
		for (int i = 0; i < closest[query].size(); i++) {
			
			string aligned = closest[query][i].seq->getAligned();
			int pos = aligned.length();
			
			//find first spot in this seq
			for (int j = aligned.length()-1; j >= 0; j--) {
				if (isalpha(aligned[j])) {
					pos = j;
					break;
				}
			}
			
			//save this spot if it is the farthest
			if (pos < rearPos) { rearPos = pos; }
		}
		
		//find last position all querySeqs[query] that is a non gap character
		aligned = querySeqs[query]->getAligned();
		pos = aligned.length();
		
		//find first spot in this seq
		for (int j = aligned.length()-1; j >= 0; j--) {
			if (isalpha(aligned[j])) {
				pos = j;
				break;
			}
		}
		
		//save this spot if it is the farthest
		if (pos < rearPos) { rearPos = pos; }

		
		//check to make sure that is not whole seq
		if ((rearPos - frontPos - 1) <= 0) {  mothurOut("Error, when I trim your sequences, the entire sequence is trimmed."); mothurOutEndLine(); exit(1);  }
		
		map<int, int> tempTrim;
		tempTrim[frontPos] = rearPos;
		
		//save trimmed locations
		trim[query] = tempTrim;
						
		//update spotMask
		map<int, int> newMap;
		int spot = 0;
		
		for (int i = frontPos; i < rearPos; i++) {
			//add to newMap
//cout << query << '\t' << i << '\t' << spotMap[query][i] << endl;
			newMap[spot] = spotMap[query][i];
			spot++;  
		}
		
		//for (it = newMap.begin(); it!=newMap.end(); it++) {
			//cout << query << '\t' << it->first << '\t' << it->second << endl;
		//}
		
		spotMap[query] = newMap;

		
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "trimSequences");
		exit(1);
	}

}
/***************************************************************************************************************/
vector<int> Ccode::findWindows(int query) {
	try {
		
		vector<int> win; 
		it = trim[query].begin();
		
		int length = it->second - it->first;
		
		//default is wanted = 10% of total length
		if (windowSizes[query] > length) { 
			mothurOut("You have slected a window larger than your sequence length after all filters, masks and trims have been done. I will use the default 10% of sequence length.");
			windowSizes[query] = length / 10;
		}else if (windowSizes[query] == 0) { windowSizes[query] = length / 10;  }
		else if (windowSizes[query] > (length / 20)) {
			mothurOut("You have selected a window that is larger than 20% of your sequence length.  This is not recommended, but I will continue anyway."); mothurOutEndLine();
		}else if (windowSizes[query] < (length / 5)) {
			mothurOut("You have selected a window that is smaller than 5% of your sequence length.  This is not recommended, but I will continue anyway."); mothurOutEndLine();
		}
		
		//save starting points of each window
		for (int m = it->first;  m < (it->second-windowSizes[query]); m+=windowSizes[query]) {  win.push_back(m);  }
		
		//save last window
		if (win[win.size()-1] < (it->first+length)) {
			win.push_back(win[win.size()-1]+windowSizes[query]); // ex. string length is 115, window is 25, without this you would get 0, 25, 50, 75
		}																									//with this you would get 1,25,50,75,100
		

		return win;
	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "findWindows");
		exit(1);
	}
}
//***************************************************************************************************************
int Ccode::getDiff(string seqA, string seqB) {
	try {
		
		int numDiff = 0;
		
		for (int i = 0; i < seqA.length(); i++) {
			//if you are both not gaps
			//if (isalpha(seqA[i]) && isalpha(seqA[i])) {
				//are you different
				if (seqA[i] != seqB[i]) { 
					 int ok; /* ok=1 means equivalent base. Checks for degenerate bases */

					/* the char in base_a and base_b have been checked and they are different */
					if ((seqA[i] == 'N') && (seqB[i] != '-')) ok = 1;
					else if ((seqB[i] == 'N') && (seqA[i] != '-')) ok = 1;
					else if ((seqA[i] == 'Y') && ((seqB[i] == 'C') || (seqB[i] == 'T'))) ok = 1;
					else if ((seqB[i] == 'Y') && ((seqA[i] == 'C') || (seqA[i] == 'T'))) ok = 1;
					else if ((seqA[i] == 'R') && ((seqB[i] == 'G') || (seqB[i] == 'A'))) ok = 1;
					else if ((seqB[i] == 'R') && ((seqA[i] == 'G') || (seqA[i] == 'A'))) ok = 1;
					else if ((seqA[i] == 'S') && ((seqB[i] == 'C') || (seqB[i] == 'G'))) ok = 1;
					else if ((seqB[i] == 'S') && ((seqA[i] == 'C') || (seqA[i] == 'G'))) ok = 1;
					else if ((seqA[i] == 'W') && ((seqB[i] == 'T') || (seqB[i] == 'A'))) ok = 1;
					else if ((seqB[i] == 'W') && ((seqA[i] == 'T') || (seqA[i] == 'A'))) ok = 1;
					else if ((seqA[i] == 'M') && ((seqB[i] == 'A') || (seqB[i] == 'C'))) ok = 1;
					else if ((seqB[i] == 'M') && ((seqA[i] == 'A') || (seqA[i] == 'C'))) ok = 1;
					else if ((seqA[i] == 'K') && ((seqB[i] == 'T') || (seqB[i] == 'G'))) ok = 1;
					else if ((seqB[i] == 'K') && ((seqA[i] == 'T') || (seqA[i] == 'G'))) ok = 1;
					else if ((seqA[i] == 'V') && ((seqB[i] == 'C') || (seqB[i] == 'A') || (seqB[i] == 'G'))) ok = 1;
					else if ((seqB[i] == 'V') && ((seqA[i] == 'C') || (seqA[i] == 'A') || (seqA[i] == 'G'))) ok = 1;
					else if ((seqA[i] == 'H') && ((seqB[i] == 'T') || (seqB[i] == 'A') || (seqB[i] == 'C'))) ok = 1;
					else if ((seqB[i] == 'H') && ((seqA[i] == 'T') || (seqA[i] == 'A') || (seqA[i] == 'C'))) ok = 1;
					else if ((seqA[i] == 'D') && ((seqB[i] == 'T') || (seqB[i] == 'A') || (seqB[i] == 'G'))) ok = 1;
					else if ((seqB[i] == 'D') && ((seqA[i] == 'T') || (seqA[i] == 'A') || (seqA[i] == 'G'))) ok = 1;
					else if ((seqA[i] == 'B') && ((seqB[i] == 'C') || (seqB[i] == 'T') || (seqB[i] == 'G'))) ok = 1;
					else if ((seqB[i] == 'B') && ((seqA[i] == 'C') || (seqA[i] == 'T') || (seqA[i] == 'G'))) ok = 1;
					else ok = 0;  /* the bases are different and not equivalent */
					
					//check if they are both blanks
					if ((seqA[i] == '.') && (seqB[i] == '-')) ok = 1;
					else if ((seqB[i] == '.') && (seqA[i] == '-')) ok = 1;
					
					if (ok == 0) {  numDiff++;  }
				}
			//}
		}
		
		return numDiff;
	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getDiff");
		exit(1);
	}
}
//***************************************************************************************************************
//tried to make this look most like ccode original implementation
void Ccode::removeBadReferenceSeqs(vector<SeqDist>& seqs, int query) {
	try {
		
		vector< vector<int> > numDiffBases;
		numDiffBases.resize(seqs.size());
		//initialize to 0
		for (int i = 0; i < numDiffBases.size(); i++) { numDiffBases[i].resize(seqs.size(),0); }
		
		it = trim[query].begin();
		int length = it->second - it->first;
		
		//calc differences from each sequence to everyother seq in the set
		for (int i = 0; i < seqs.size(); i++) {
			
			string seqA = seqs[i].seq->getAligned().substr(it->first, length);
			
			//so you don't calc i to j and j to i since they are the same
			for (int j = 0; j < i; j++) {
				
				string seqB = seqs[j].seq->getAligned().substr(it->first, length);
				
				//compare strings
				int numDiff = getDiff(seqA, seqB);
				
				numDiffBases[i][j] = numDiff;
				numDiffBases[j][i] = numDiff;
			}
		}
		
		//initailize remove to 0
		vector<int> remove;  remove.resize(seqs.size(), 0);
		float top = ((20*length) / (float) 100);
		float bottom = ((0.5*length) / (float) 100);
		
		//check each numDiffBases and if any are higher than threshold set remove to 1 so you can remove those seqs from the closest set
		for (int i = 0; i < numDiffBases.size(); i++) {
			for (int j = 0; j < i; j++) {	
				//are you more than 20% different
				if (numDiffBases[i][j] > top)		{  remove[j] = 1;  }
				//are you less than 0.5% different
				if (numDiffBases[i][j] < bottom)	{  remove[j] = 1;  }
			}
		}
		
		int numSeqsLeft = 0;
		
		//count seqs that are not going to be removed
		for (int i = 0; i < remove.size(); i++) {  
			if (remove[i] == 0)  { numSeqsLeft++;  }
		}
		
		//if you have enough then remove bad ones
		if (numSeqsLeft >= 3) {
			vector<SeqDist> goodSeqs;
			//remove bad seqs
			for (int i = 0; i < remove.size(); i++) {
				if (remove[i] == 0) { 
					goodSeqs.push_back(seqs[i]);
				}
			}
			
			seqs = goodSeqs;
			
		}else { //warn, but dont remove any
			mothurOut(querySeqs[query]->getName() + " does not have an adaquate number of reference sequences that are within 20% and 0.5% similarity.  I will continue, but please check."); mothurOutEndLine();  
		}
			

	}
	catch(exception& e) {
		errorOut(e, "Ccode", "removeBadReferenceSeqs");
		exit(1);
	}
}
//***************************************************************************************************************
vector< vector<SeqDist> > Ccode::findClosest(int start, int end, int numWanted) {
	try{
	
		vector< vector<SeqDist> > topMatches;  topMatches.resize(querySeqs.size());
	
		float smallestOverall, smallestLeft, smallestRight;
		smallestOverall = 1000;  smallestLeft = 1000;  smallestRight = 1000;
		
		//for each sequence in querySeqs - find top matches to use as reference
		for(int j = start; j < end; j++){
			
			Sequence query = *(querySeqs[j]);
			
			vector<SeqDist> distances;
			
			//calc distance to each sequence in template seqs
			for (int i = 0; i < templateSeqs.size(); i++) {
			
				Sequence ref = *(templateSeqs[i]); 
					
				//find overall dist
				distCalc->calcDist(query, ref);
				float dist = distCalc->getDist();	
				
				//save distance
				SeqDist temp;
				temp.seq = templateSeqs[i];
				temp.dist = dist;

				distances.push_back(temp);
			}
			
			sort(distances.begin(), distances.end(), compareSeqDist);
			
			//save the number of top matches wanted
			for (int h = 0; h < numWanted; h++) {
				topMatches[j].push_back(distances[h]);
			}
		}
			
		return topMatches;

	}
	catch(exception& e) {
		errorOut(e, "Ccode", "findClosestSides");
		exit(1);
	}
}
/**************************************************************************************************/
//find the distances from each reference sequence to every other reference sequence for each window for this query
void Ccode::getAverageRef(vector<SeqDist> ref, int query) {
	try {
		
		vector< vector< vector<int> > > diffs;  //diffs[0][1][2] is the number of differences between ref seq 0 and ref seq 1 at window 2.
		
		//initialize diffs vector
		diffs.resize(ref.size());
		for (int i = 0; i < diffs.size(); i++) {  
			diffs[i].resize(ref.size());
			for (int j = 0; j < diffs[i].size(); j++) {
				diffs[i][j].resize(windows[query].size(), 0);
			}
		}
		
		it = trim[query].begin();
				
		//find the distances from each reference sequence to every other reference sequence for each window for this query		
		for (int i = 0; i < ref.size(); i++) {
			
			string refI = ref[i].seq->getAligned();
			
			//j<i, so you don't find distances from i to j and then j to i.
			for (int j = 0; j < i; j++) {
			
				string refJ = ref[j].seq->getAligned();
			
				for (int k = 0; k < windows[query].size(); k++) {
					
					string refIWindowk, refJWindowk;
					
					if (k < windows[query].size()-1) {
						//get window strings
						refIWindowk = refI.substr(windows[query][k], windowSizes[query]);
						refJWindowk = refJ.substr(windows[query][k], windowSizes[query]);
					}else { //last window may be smaller than rest - see findwindows
						//get window strings
						refIWindowk = refI.substr(windows[query][k], (it->second-windows[query][k]));
						refJWindowk = refJ.substr(windows[query][k], (it->second-windows[query][k]));
					}
					
					//find differences
					int diff = getDiff(refIWindowk, refJWindowk);
//cout <<  i << '\t' << j << '\t' << k << '\t' << diff << endl;
					//save differences in [i][j][k] and [j][i][k] since they are the same
					diffs[i][j][k] = diff;
					diffs[j][i][k] = diff;

				}//k
				
			}//j
		
		}//i
		
		//initialize sumRef for this query 
		sumRef[query].resize(windows[query].size(), 0);
		sumSquaredRef[query].resize(windows[query].size(), 0);
		averageRef[query].resize(windows[query].size(), 0);
		
		//find the sum of the differences for hte reference sequences
		for (int i = 0; i < diffs.size(); i++) {
			for (int j = 0; j < i; j++) {
			
				//increment this querys reference sequences combos
				refCombo[query]++;
				
				for (int k = 0; k < diffs[i][j].size(); k++) {
					sumRef[query][k] += diffs[i][j][k];
					sumSquaredRef[query][k] += (diffs[i][j][k]*diffs[i][j][k]);
				}//k
				
			}//j
		}//i

		
		//find the average of the differences for the references for each window
		for (int i = 0; i < windows[query].size(); i++) {
			averageRef[query][i] = sumRef[query][i] / (float) refCombo[query];
		}
		
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getAverageRef");
		exit(1);
	}
}
/**************************************************************************************************/
void Ccode::getAverageQuery (vector<SeqDist> ref, int query) {
	try {
	
		vector< vector<int> >  diffs;  //diffs[1][2] is the number of differences between querySeqs[query] and ref seq 1 at window 2.
		
		//initialize diffs vector
		diffs.resize(ref.size());
		for (int j = 0; j < diffs.size(); j++) {
			diffs[j].resize(windows[query].size(), 0);
		}
		
		it = trim[query].begin();
							
		string refQuery = querySeqs[query]->getAligned();
			 
		//j<i, so you don't find distances from i to j and then j to i.
		for (int j = 0; j < ref.size(); j++) {
			 
			 string refJ = ref[j].seq->getAligned();
			 
			 for (int k = 0; k < windows[query].size(); k++) {
					
					string QueryWindowk, refJWindowk;
					
					if (k < windows[query].size()-1) {
						//get window strings
						QueryWindowk = refQuery.substr(windows[query][k], windowSizes[query]);
						refJWindowk = refJ.substr(windows[query][k], windowSizes[query]);					
					}else { //last window may be smaller than rest - see findwindows
						//get window strings
						QueryWindowk = refQuery.substr(windows[query][k], (it->second-windows[query][k]));
						refJWindowk = refJ.substr(windows[query][k], (it->second-windows[query][k]));
					}
					
					//find differences
					int diff = getDiff(QueryWindowk, refJWindowk);
//cout  << j << '\t' << k << '\t' << diff << endl;						
					//save differences 
					diffs[j][k] = diff;
			 
			 }//k
		}//j
			 
		
		//initialize sumRef for this query 
		sumQuery[query].resize(windows[query].size(), 0);
		sumSquaredQuery[query].resize(windows[query].size(), 0);
		averageQuery[query].resize(windows[query].size(), 0);
		
		//find the sum of the differences 
		for (int j = 0; j < diffs.size(); j++) {
			for (int k = 0; k < diffs[j].size(); k++) {
				sumQuery[query][k] += diffs[j][k];
				sumSquaredQuery[query][k] += (diffs[j][k]*diffs[j][k]);
			}//k
		}//j
		
		
		//find the average of the differences for the references for each window
		for (int i = 0; i < windows[query].size(); i++) {
			averageQuery[query][i] = sumQuery[query][i] / (float) ref.size();
		}

	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getAverageQuery");
		exit(1);
	}
}
/**************************************************************************************************/
void Ccode::findVarianceRef (int query) {
	try {
		
		varRef[query].resize(windows[query].size(), 0);
		sdRef[query].resize(windows[query].size(), 0);
		
		//for each window
		for (int i = 0; i < windows[query].size(); i++) {
			varRef[query][i] = (sumSquaredRef[query][i] - ((sumRef[query][i]*sumRef[query][i])/(float)refCombo[query])) / (float)(refCombo[query]-1);
			sdRef[query][i] = sqrt(varRef[query][i]);
			
			//set minimum error rate to 0.001 - to avoid potential divide by zero - not sure if this is necessary but it follows ccode implementation
			if (averageRef[query][i] < 0.001)			{	averageRef[query][i] = 0.001;		}
			if (sumRef[query][i] < 0.001)				{	sumRef[query][i] = 0.001;			}
			if (varRef[query][i] < 0.001)				{	varRef[query][i] = 0.001;			}
			if (sumSquaredRef[query][i] < 0.001)		{	sumSquaredRef[query][i] = 0.001;	}
			if (sdRef[query][i] < 0.001)				{	sdRef[query][i] = 0.001;			}
				
		}
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "findVarianceRef");
		exit(1);
	}
}
/**************************************************************************************************/
void Ccode::findVarianceQuery (int query) {
	try {
		varQuery[query].resize(windows[query].size(), 0);
		sdQuery[query].resize(windows[query].size(), 0);
		
		//for each window
		for (int i = 0; i < windows[query].size(); i++) {
			varQuery[query][i] = (sumSquaredQuery[query][i] - ((sumQuery[query][i]*sumQuery[query][i])/(float) closest[query].size())) / (float) (closest[query].size()-1);
			sdQuery[query][i] = sqrt(varQuery[query][i]);
			
			//set minimum error rate to 0.001 - to avoid potential divide by zero - not sure if this is necessary but it follows ccode implementation
			if (averageQuery[query][i] < 0.001)			{	averageQuery[query][i] = 0.001;		}
			if (sumQuery[query][i] < 0.001)				{	sumQuery[query][i] = 0.001;			}
			if (varQuery[query][i] < 0.001)				{	varQuery[query][i] = 0.001;			}
			if (sumSquaredQuery[query][i] < 0.001)		{	sumSquaredQuery[query][i] = 0.001;	}
			if (sdQuery[query][i] < 0.001)				{	sdQuery[query][i] = 0.001;			}
		}

	}
	catch(exception& e) {
		errorOut(e, "Ccode", "findVarianceQuery");
		exit(1);
	}
}
/**************************************************************************************************/
void Ccode::determineChimeras (int query) {
	try {
		
		isChimericConfidence[query].resize(windows[query].size(), false);
		isChimericTStudent[query].resize(windows[query].size(), false);
		isChimericANOVA[query].resize(windows[query].size(), false);
		anova[query].resize(windows[query].size());

		
		//for each window
		for (int i = 0; i < windows[query].size(); i++) {
		
			//get confidence limits
			float t = getT(closest[query].size()-1);  //how many seqs you are comparing to this querySeq
			float dsUpper = (averageQuery[query][i] + (t * sdQuery[query][i])) / averageRef[query][i]; 
			float dsLower = (averageQuery[query][i] - (t * sdQuery[query][i])) / averageRef[query][i]; 
//cout << t << '\t' << "ds upper = " << dsUpper << " dsLower = " << dsLower << endl;			
			if ((dsUpper > 1.0) && (dsLower > 1.0) && (averageQuery[query][i] > averageRef[query][i])) {  /* range does not include 1 */
					isChimericConfidence[query][i] = true;   /* significantly higher at P<0.05 */ 
//cout << i << " made it here" << endl;
			}
			
			//student t test
			int degreeOfFreedom = refCombo[query] + closest[query].size() - 2;
			float denomForT = (((refCombo[query]-1) * varQuery[query][i] + (closest[query].size() - 1) * varRef[query][i]) / (float) degreeOfFreedom) * ((refCombo[query] + closest[query].size()) / (float) (refCombo[query] * closest[query].size()));	/* denominator, without sqrt(), for ts calculations */
				
			float ts = fabs((averageQuery[query][i] - averageRef[query][i]) / (sqrt(denomForT)));  /* value of ts for t-student test */	
			t = getT(degreeOfFreedom);
//cout << i << '\t' << t << '\t' << ts << endl;			
			if ((ts >= t) && (averageQuery[query][i] > averageRef[query][i])) {   
				isChimericTStudent[query][i] = true;   /* significantly higher at P<0.05 */ 
			}
		
			//anova test
			float value1 = sumQuery[query][i] + sumRef[query][i];
			float value2 = sumSquaredQuery[query][i] + sumSquaredRef[query][i];
			float value3 = ((sumQuery[query][i]*sumQuery[query][i]) / (float) (closest[query].size())) + ((sumRef[query][i] * sumRef[query][i]) / (float) refCombo[query]);
			float value4 = (value1 * value1) / ( (float) (closest[query].size() + refCombo[query]) );
			float value5 = value2 - value4;
			float value6 = value3 - value4;
			float value7 = value5 - value6;
			float value8 = value7 / ((float) degreeOfFreedom);
			float anovaValue = value6 / value8;
			
			float f = getF(degreeOfFreedom);
			
			 if ((anovaValue >= f) && (averageQuery[query][i] > averageRef[query][i]))  {
	               isChimericANOVA[query][i] = true;   /* significant P<0.05 */
	        }
			
			if (isnan(anovaValue) || isinf(anovaValue)) { anovaValue = 0.0; }
			
			anova[query][i] = anovaValue;
		}
		
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "determineChimeras");
		exit(1);
	}
}
/**************************************************************************************************/	
float Ccode::getT(int numseq) {
	try {
	
		float tvalue = 0;
		
		/* t-student critical values for different degrees of freedom and alpha 0.1 in one-tail tests (equivalent to 0.05) */
		if (numseq > 120) tvalue = 1.645;
		else if (numseq > 60) tvalue = 1.658;
        else if (numseq > 40) tvalue = 1.671;
        else if (numseq > 30) tvalue = 1.684;
        else if (numseq > 29) tvalue = 1.697;
        else if (numseq > 28) tvalue = 1.699;
        else if (numseq > 27) tvalue = 1.701;
        else if (numseq > 26) tvalue = 1.703;
        else if (numseq > 25) tvalue = 1.706;
        else if (numseq > 24) tvalue = 1.708;
        else if (numseq > 23) tvalue = 1.711;
        else if (numseq > 22) tvalue = 1.714;
        else if (numseq > 21) tvalue = 1.717;
        else if (numseq > 20) tvalue = 1.721;
        else if (numseq > 19) tvalue = 1.725;
        else if (numseq > 18) tvalue = 1.729;
        else if (numseq > 17) tvalue = 1.734;
        else if (numseq > 16) tvalue = 1.740;
        else if (numseq > 15) tvalue = 1.746;
        else if (numseq > 14) tvalue = 1.753;
        else if (numseq > 13) tvalue = 1.761;
        else if (numseq > 12) tvalue = 1.771;
        else if (numseq > 11) tvalue = 1.782;
        else if (numseq > 10) tvalue = 1.796;
        else if (numseq > 9) tvalue = 1.812;
        else if (numseq > 8) tvalue = 1.833;
        else if (numseq > 7) tvalue = 1.860;
        else if (numseq > 6) tvalue = 1.895;
        else if (numseq > 5) tvalue = 1.943;
        else if (numseq > 4) tvalue = 2.015;
        else if (numseq > 3) tvalue = 2.132;
        else if (numseq > 2) tvalue = 2.353;
        else if (numseq > 1) tvalue = 2.920;
		else if (numseq <= 1) {
			mothurOut("Two or more reference sequences are required, your data will be flawed.\n"); mothurOutEndLine();
		}
		
		return tvalue;
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getT");
		exit(1);
	}
}
/**************************************************************************************************/	
float Ccode::getF(int numseq) {
	try {
	
		float fvalue = 0;
		
		 /* F-Snedecor critical values for v1=1 and different degrees of freedom v2 and alpha 0.05 */
        if (numseq > 120) fvalue = 3.84;
        else if (numseq > 60) fvalue = 3.92;
        else if (numseq > 40) fvalue = 4.00;
        else if (numseq > 30) fvalue = 4.08;
        else if (numseq > 29) fvalue = 4.17;
        else if (numseq > 28) fvalue = 4.18;
        else if (numseq > 27) fvalue = 4.20;
        else if (numseq > 26) fvalue = 4.21;
        else if (numseq > 25) fvalue = 4.23;
        else if (numseq > 24) fvalue = 4.24;
        else if (numseq > 23) fvalue = 4.26;
        else if (numseq > 22) fvalue = 4.28;
        else if (numseq > 21) fvalue = 4.30;
        else if (numseq > 20) fvalue = 4.32;
        else if (numseq > 19) fvalue = 4.35;
        else if (numseq > 18) fvalue = 4.38;
        else if (numseq > 17) fvalue = 4.41;
        else if (numseq > 16) fvalue = 4.45;
        else if (numseq > 15) fvalue = 4.49;
        else if (numseq > 14) fvalue = 4.54;
        else if (numseq > 13) fvalue = 4.60;
        else if (numseq > 12) fvalue = 4.67;
        else if (numseq > 11) fvalue = 4.75;
        else if (numseq > 10) fvalue = 4.84;
        else if (numseq > 9) fvalue = 4.96;
        else if (numseq > 8) fvalue = 5.12;
        else if (numseq > 7) fvalue = 5.32;
        else if (numseq > 6) fvalue = 5.59;
        else if (numseq > 5) fvalue = 5.99;
        else if (numseq > 4) fvalue = 6.61;
        else if (numseq > 3) fvalue = 7.71;
        else if (numseq > 2) fvalue = 10.1;
        else if (numseq > 1) fvalue = 18.5;
        else if (numseq > 0) fvalue = 161;
		else if (numseq <= 0) {
			mothurOut("Two or more reference sequences are required, your data will be flawed.\n"); mothurOutEndLine();
        }
		
		return fvalue;
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getF");
		exit(1);
	}
}

/**************************************************************************************************/
void Ccode::createProcessesClosest() {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		vector<int> processIDS;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
				
				mothurOut("Finding top matches for sequences " + toString(lines[process]->start) + " to " + toString(lines[process]->end)); mothurOutEndLine();
				closest = findClosest(lines[process]->start, lines[process]->end, numWanted);
				mothurOut("Done finding top matches for sequences " +  toString(lines[process]->start) + " to " + toString(lines[process]->end)); mothurOutEndLine();
				
				//write out data to file so parent can read it
				ofstream out;
				string s = toString(getpid()) + ".temp";
				openOutputFile(s, out);
							
				//output pairs
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					 for (int j = 0; j < closest[i].size(); j++) {
						closest[i][j].seq->printSequence(out);
					 }
				}
				out << ">" << endl; //to stop sequence read
				
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					 for (int j = 0; j < closest[i].size(); j++) {
						out << closest[i][j].dist << '\t';
					 }
					 out << endl;
				}
				out.close();
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//get data created by processes
		for (int i=0;i<processors;i++) { 
			ifstream in;
			string s = toString(processIDS[i]) + ".temp";
			openInputFile(s, in);
			
			vector< vector<Sequence*> > tempClosest;  tempClosest.resize(querySeqs.size());
			//get pairs
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				vector<Sequence*> tempVector;
	
				for (int j = 0; j < numWanted; j++) {
				
					Sequence* temp = new Sequence(in);
					gobble(in);
						
					tempVector.push_back(temp);
				}
				
				tempClosest[k] = tempVector;
			}
			
			string junk;
			in >> junk;  gobble(in);  // to get ">"
		
			vector< vector<float> > dists; dists.resize(querySeqs.size());
			
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				dists[k].resize(numWanted);
				for (int j = 0; j < numWanted; j++) {
					in >> dists[k][j];
				}
				gobble(in);
			
			} 

			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				closest[k].resize(numWanted);
				for (int j = 0; j < closest[k].size(); j++) {
					closest[k][j].seq = tempClosest[k][j];
					closest[k][j].dist = dists[k][j]; 
				}
			} 

			in.close();
			remove(s.c_str());
		}
			
	
#else
		closest = findClosest(lines[0]->start, lines[0]->end, numWanted);
#endif	

	}
	catch(exception& e) {
		errorOut(e, "Ccode", "createProcessesClosest");
		exit(1);
	}
}

//***************************************************************************************************************
void Ccode::createProcessesRemoveBad() {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		vector<int> processIDS;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
							
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					removeBadReferenceSeqs(closest[i], i);
				}
				
				//write out data to file so parent can read it
				ofstream out;
				string s = toString(getpid()) + ".temp";
				openOutputFile(s, out);
				
				//output pairs
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					out << closest[i].size() << endl;
					 for (int j = 0; j < closest[i].size(); j++) {
						closest[i][j].seq->printSequence(out);
					 }
					 out << ">" << endl; //to stop sequence read
				}
				
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					 for (int j = 0; j < closest[i].size(); j++) {
						out << closest[i][j].dist << '\t';
					 }
					 out << endl;
				}

				out.close();
				
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
	
		//get data created by processes
		for (int i=0;i<processors;i++) { 
			ifstream in;
			string s = toString(processIDS[i]) + ".temp";
			openInputFile(s, in);
		
			vector< vector<Sequence*> > tempClosest;  tempClosest.resize(querySeqs.size());
			vector<int> sizes; 
			//get pairs
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
		
				int num;
				in >> num; 
				sizes.push_back(num);
				gobble(in);
				
				vector<Sequence*> tempVector;
				
				for (int j = 0; j < num; j++) {
				
					Sequence* temp = new Sequence(in);
					gobble(in);
						
					tempVector.push_back(temp);
				}
				string junk;
				in >> junk;  gobble(in);  // to get ">"

				tempClosest[k] = tempVector;
			}
			
			vector< vector<float> > dists; dists.resize(querySeqs.size());
			int count = 0;
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				dists[k].resize(sizes[count]);
				for (int j = 0; j < sizes[count]; j++) {
					in >> dists[k][j];
				}
				gobble(in);
				count++;
			} 
			
			count = 0;
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				for (int j = 0; j < sizes[count]; j++) {
					closest[k][j].seq = tempClosest[k][j];
					closest[k][j].dist = dists[k][j]; 
				}
				count++;
			} 

			in.close();
			remove(s.c_str());
		}
#else
		for (int i = 0; i < closest.size(); i++) {
			removeBadReferenceSeqs(closest[i], i);
		}
#endif		
	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "createProcessesRemoveBad");
		exit(1);
	}
}
//***************************************************************************************************************
void Ccode::createProcessesAverages() {
	try {
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		vector<int> processIDS;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
							
				//find the averages for each querys references 
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					getAverageRef(closest[i], i);  //fills sumRef[i], averageRef[i], sumSquaredRef[i] and refCombo[i].
				}
				
				//find the averages for the query 
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					getAverageQuery(closest[i], i);  //fills sumQuery[i], averageQuery[i], sumSquaredQuery[i].
				}
				
				
				//write out data to file so parent can read it
				ofstream out;
				string s = toString(getpid()) + ".temp";
				openOutputFile(s, out);
				
				//output pairs
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					for (int j = 0; j < windows[i].size(); j++) {
						out << sumRef[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << averageRef[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << sumSquaredRef[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << sumQuery[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << averageQuery[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << sumSquaredQuery[i][j] << '\t';
					}
					out << endl;
					out << refCombo[i] << endl;
				}

				out.close();
				
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
	
		//get data created by processes
		for (int i=0;i<processors;i++) { 
			ifstream in;
			string s = toString(processIDS[i]) + ".temp";
			openInputFile(s, in);
		
			//get pairs
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				sumRef[k].resize(windows[k].size());
				averageRef[k].resize(windows[k].size());
				sumSquaredRef[k].resize(windows[k].size());
				averageQuery[k].resize(windows[k].size());
				sumQuery[k].resize(windows[k].size());
				sumSquaredQuery[k].resize(windows[k].size());
				
				for (int j = 0; j < windows[k].size(); j++) {
					in >> sumRef[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> averageRef[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> sumSquaredRef[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> sumQuery[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> averageQuery[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> sumSquaredQuery[k][j];
				}
				gobble(in);
				in >> refCombo[k];
				gobble(in);
			}
			
			in.close();
			remove(s.c_str());
		}
#else
		//find the averages for each querys references 
		for (int i = 0; i < querySeqs.size(); i++) {
			getAverageRef(closest[i], i);  //fills sumRef[i], averageRef[i], sumSquaredRef[i] and refCombo[i].
		}
	
		//find the averages for the query 
		for (int i = 0; i < querySeqs.size(); i++) {
			getAverageQuery(closest[i], i);  //fills sumQuery[i], averageQuery[i], sumSquaredQuery[i].
		}

#endif		
	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "createProcessesAverages");
		exit(1);
	}
}
//***************************************************************************************************************
void Ccode::createProcessesVariances() {
	try {
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		vector<int> processIDS;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
							
				//find the averages for each querys references 
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					findVarianceRef(i);  //fills varRef[i] and sdRef[i] also sets minimum error rate to 0.001 to avoid divide by 0.
				}
				
				//find the averages for the query 
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					findVarianceQuery(i);  //fills varQuery[i] and sdQuery[i] also sets minimum error rate to 0.001 to avoid divide by 0.
				}
				
				
				//write out data to file so parent can read it
				ofstream out;
				string s = toString(getpid()) + ".temp";
				openOutputFile(s, out);
				
				//output pairs
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					for (int j = 0; j < windows[i].size(); j++) {
						out << varRef[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << sdRef[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << varQuery[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << sdQuery[i][j] << '\t';
					}
					out << endl;
				}

				out.close();
				
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
	
		//get data created by processes
		for (int i=0;i<processors;i++) { 
			ifstream in;
			string s = toString(processIDS[i]) + ".temp";
			openInputFile(s, in);
		
			//get pairs
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				varRef[k].resize(windows[k].size());
				sdRef[k].resize(windows[k].size());
				varQuery[k].resize(windows[k].size());
				sdQuery[k].resize(windows[k].size());
								
				for (int j = 0; j < windows[k].size(); j++) {
					in >> varRef[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> sdRef[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> varQuery[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> sdQuery[k][j];
				}
				gobble(in);
			}
			
			in.close();
			remove(s.c_str());
		}
#else
			//find the averages for each querys references 
			for (int i = 0; i < querySeqs.size(); i++) {
				findVarianceRef(i);  //fills varRef[i] and sdRef[i] also sets minimum error rate to 0.001 to avoid divide by 0.
			}
			
			//find the averages for the query 
			for (int i = 0; i < querySeqs.size(); i++) {
				findVarianceQuery(i);  //fills varQuery[i] and sdQuery[i] also sets minimum error rate to 0.001 to avoid divide by 0.
			}
#endif		
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "createProcessesVariances");
		exit(1);
	}
}
//***************************************************************************************************************
void Ccode::createProcessesDetermine() {
	try {
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		vector<int> processIDS;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
							
				//find the averages for each querys references 
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					determineChimeras(i);  //fills anova, isChimericConfidence[i], isChimericTStudent[i] and isChimericANOVA[i]. 
				}

				//write out data to file so parent can read it
				ofstream out;
				string s = toString(getpid()) + ".temp";
				openOutputFile(s, out);
				
				//output pairs
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					for (int j = 0; j < windows[i].size(); j++) {
						out << anova[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << isChimericConfidence[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << isChimericTStudent[i][j] << '\t';
					}
					out << endl;
					for (int j = 0; j < windows[i].size(); j++) {
						out << isChimericANOVA[i][j] << '\t';
					}
					out << endl;
				}

				out.close();
				
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
	
		//get data created by processes
		for (int i=0;i<processors;i++) { 
			ifstream in;
			string s = toString(processIDS[i]) + ".temp";
			openInputFile(s, in);
		
			//get pairs
			for (int k = lines[i]->start; k < lines[i]->end; k++) {
				anova[k].resize(windows[k].size());
				isChimericConfidence[k].resize(windows[k].size(), false);
				isChimericTStudent[k].resize(windows[k].size(), false);
				isChimericANOVA[k].resize(windows[k].size(), false);
				int tempBool;
								
				for (int j = 0; j < windows[k].size(); j++) {
					in >> anova[k][j];
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> tempBool;
					if (tempBool == 1) {  isChimericConfidence[k][j] = true;  }
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> tempBool;
					if (tempBool == 1) {  isChimericTStudent[k][j] = true;  }
				}
				gobble(in);
				for (int j = 0; j < windows[k].size(); j++) {
					in >> tempBool;
					if (tempBool == 1) {  isChimericANOVA[k][j] = true;  }
				}
				gobble(in);
			}
			
			in.close();
			remove(s.c_str());
		}
	#else
			for (int i = 0; i < querySeqs.size(); i++) {
				determineChimeras(i);  //fills anova, isChimericConfidence[i], isChimericTStudent[i] and isChimericANOVA[i].
			}
	#endif		

	}
	catch(exception& e) {
		errorOut(e, "Ccode", "createProcessesDetermine");
		exit(1);
	}
}
//***************************************************************************************************************



