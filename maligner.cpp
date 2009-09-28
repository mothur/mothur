/*
 *  maligner.cpp
 *  Mothur
 *
 *  Created by westcott on 9/23/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "maligner.h"

/***********************************************************************/
Maligner::Maligner(vector<Sequence*> temp, int num, int match, int misMatch, float div, int minCov) :
		db(temp), numWanted(num), matchScore(match), misMatchPenalty(misMatch), minDivR(div), minCoverage(minCov) {}
/***********************************************************************/
string Maligner::getResults(Sequence* q) {
	try {
		
		//make copy so trimming doesn't destroy query from calling class - remember to deallocate
		query = new Sequence(q->getName(), q->getAligned());
		
		string chimera;
	
		decalc = new DeCalculator();
		
		//find closest seqs to query in template - returns copies of seqs so trim does not destroy - remember to deallocate
		refSeqs = decalc->findClosest(query, db, numWanted);
		
		ofstream out;
		string outFile = "parentsOf" + query->getName();
		openOutputFile(outFile, out);
		for (int i = 0; i < refSeqs.size(); i++) {   refSeqs[i]->printSequence(out);	}
		out.close();
		
		refSeqs = minCoverageFilter(refSeqs);
		
		if (refSeqs.size() < 2)  { 
			for (int i = 0; i < refSeqs.size(); i++) {  delete refSeqs[i];	}
			percentIdenticalQueryChimera = 0.0;
			return "unknown"; 
		}
		
		int chimeraPenalty = computeChimeraPenalty();
	
		//trims seqs to first non gap char in all seqs and last non gap char in all seqs
		decalc->trimSeqs(query, refSeqs);
		
		vector<Sequence*> temp = refSeqs;
		temp.push_back(query);
		
		verticalFilter(temp);

		vector< vector<score_struct> > matrix = buildScoreMatrix(query->getAligned().length(), refSeqs.size()); //builds and initializes
		
		fillScoreMatrix(matrix, refSeqs, chimeraPenalty);
		
		vector<score_struct> path = extractHighestPath(matrix);
		
		vector<trace_struct> trace = mapTraceRegionsToAlignment(path, refSeqs);
		
		if (trace.size() > 1) {		chimera = "yes";	}
		else { chimera = "no";	}
		
		int traceStart = path[0].col;
		int traceEnd = path[path.size()-1].col;
		
		string queryInRange = query->getAligned();
		queryInRange = queryInRange.substr(traceStart, (traceEnd-traceStart+1));
	
		string chimeraSeq = constructChimericSeq(trace, refSeqs);
		
		percentIdenticalQueryChimera = computePercentID(queryInRange, chimeraSeq);
		
		delete decalc;
		
		//save output results
		for (int i = 0; i < trace.size(); i++) {
			int regionStart = trace[i].col;
			int regionEnd = trace[i].oldCol;
			int seqIndex = trace[i].row;
			
			results temp;
			
			temp.parent = refSeqs[seqIndex]->getName();
			temp.regionStart = regionStart;
			temp.regionEnd = regionEnd;
			
			string parentInRange = refSeqs[seqIndex]->getAligned();
			parentInRange = parentInRange.substr(traceStart, (traceEnd-traceStart+1));
			
			temp.queryToParent = computePercentID(queryInRange, parentInRange);
			temp.divR = (percentIdenticalQueryChimera / temp.queryToParent);
			
			string queryInRegion = query->getAligned();
			queryInRegion = queryInRegion.substr(regionStart, (regionEnd-regionStart+1));
			
			string parentInRegion = refSeqs[seqIndex]->getAligned();
			parentInRegion = parentInRegion.substr(regionStart, (regionEnd-regionStart+1));
			
			temp.queryToParentLocal = computePercentID(queryInRegion, parentInRegion);
		
			outputResults.push_back(temp);
		}
		
		//free memory
		delete query;
		for (int i = 0; i < refSeqs.size(); i++) {  delete refSeqs[i];	}
		
		return chimera;
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "getResults");
		exit(1);
	}
}
/***********************************************************************/
//removes top matches that do not have minimum coverage with query.
vector<Sequence*> Maligner::minCoverageFilter(vector<Sequence*> ref){  
	try {
		vector<Sequence*> newRefs;
		
		string queryAligned = query->getAligned();
		
		for (int i = 0; i < ref.size(); i++) {
			
			string refAligned = ref[i]->getAligned();
			
			int numBases = 0;
			int numCovered = 0;
			
			//calculate coverage
			for (int j = 0; j < queryAligned.length(); j++) {
				
				if (isalpha(queryAligned[j])) {
					numBases++;
					
					if (isalpha(refAligned[j])) {
						numCovered++;
					}
				}
			}
			
			int coverage = ((numCovered/(float)numBases)*100);
			
			//if coverage above minimum
			if (coverage > minCoverage) {
				newRefs.push_back(ref[i]);
			}
		}
		
		return newRefs;
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "minCoverageFilter");
		exit(1);
	}
}
/***********************************************************************/
// a breakpoint should yield fewer mismatches than this number with respect to the best parent sequence.
int Maligner::computeChimeraPenalty() {
	try {
		
		int numAllowable = ((1.0 - (1.0/minDivR)) * query->getNumBases());
	
		int penalty = int(numAllowable + 1) * misMatchPenalty;
											 
		return penalty;

	}
	catch(exception& e) {
		errorOut(e, "Maligner", "computeChimeraPenalty");
		exit(1);
	}
}
/***********************************************************************/
//this is a vertical filter
void Maligner::verticalFilter(vector<Sequence*> seqs) {
	try {
		vector<int> gaps;	gaps.resize(query->getAligned().length(), 0);
		
		string filterString = (string(query->getAligned().length(), '1'));
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is a gap
				if ((seqAligned[j] == '-') || (seqAligned[j] == '.'))	{	gaps[j]++;	}
			}
		}
		
		//zero out spot where all sequences have blanks
		int numColRemoved = 0;
		for(int i = 0; i < seqs[0]->getAligned().length(); i++){
			if(gaps[i] == seqs.size())	{	filterString[i] = '0'; 	numColRemoved++;  }
		}
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			string newAligned = "";
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is not a gap
				if (filterString[j] == '1') { newAligned += seqAligned[j]; }
			}
			
			seqs[i]->setAligned(newAligned);
		}

	
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "verticalFilter");
		exit(1);
	}
}
//***************************************************************************************************************
vector< vector<score_struct> > Maligner::buildScoreMatrix(int cols, int rows) {
	try{
		
		vector< vector<score_struct> > m; m.resize(rows);
		
		for (int i = 0; i < m.size(); i++) {
			for (int j = 0; j < cols; j++) {
				
				//initialize each cell
				score_struct temp;
				temp.prev = -1;
				temp.score = -9999999;
				temp.col = j;
				temp.row = i;
				
				m[i].push_back(temp);
			}
		}
		
		return m;
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "buildScoreMatrix");
		exit(1);
	}
}
//***************************************************************************************************************
void Maligner::fillScoreMatrix(vector<vector<score_struct> >& m, vector<Sequence*> seqs, int penalty) {
	try{
		
		//get matrix dimensions
		int numCols = query->getAligned().length();
		int numRows = seqs.size();
		
		//initialize first col
		string queryAligned = query->getAligned();
		for (int i = 0; i < numRows; i++) {
			string subjectAligned = seqs[i]->getAligned();
			
			//are you both gaps?
			if ((!isalpha(queryAligned[0])) && (!isalpha(subjectAligned[0]))) {
				m[i][0].score = 0;
			}else if (queryAligned[0] == subjectAligned[0]) {
				m[i][0].score = matchScore;
			}else{
				m[i][0].score = 0;
			}
		}
		
		//fill rest of matrix
		for (int j = 1; j < numCols; j++) {  //iterate through matrix columns
		
			for (int i = 0; i < numRows; i++) {  //iterate through matrix rows
				
				string subjectAligned = seqs[i]->getAligned();
				
				int matchMisMatchScore = 0;
				//are you both gaps?
				if ((!isalpha(queryAligned[j])) && (!isalpha(subjectAligned[j]))) {
					//leave the same
				}else if ((toupper(queryAligned[j]) == 'N') || (toupper(subjectAligned[j]) == 'N')) {
					//leave the same
				}else if (queryAligned[j] == subjectAligned[j]) {
					matchMisMatchScore = matchScore;
				}else if (queryAligned[j] != subjectAligned[j]) {
					matchMisMatchScore = misMatchPenalty;
				}
				
				//compute score based on previous columns scores
				for (int prevIndex = 0; prevIndex < numRows; prevIndex++) { //iterate through rows
					
					int sumScore = matchMisMatchScore + m[prevIndex][j-1].score;
					
					//you are not at yourself
					if (prevIndex != i) {   sumScore += penalty;	}
					if (sumScore < 0)	{	sumScore = 0;			}
					
					if (sumScore > m[i][j].score) {
						m[i][j].score = sumScore;
						m[i][j].prev = prevIndex;
					}
				}
			}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "fillScoreMatrix");
		exit(1);
	}
}
//***************************************************************************************************************
vector<score_struct> Maligner::extractHighestPath(vector<vector<score_struct> > m) {
	try {
	
		//get matrix dimensions
		int numCols = query->getAligned().length();
		int numRows = m.size();
	
	
		//find highest score scoring matrix
		score_struct highestStruct;
		int highestScore = 0;
		
		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				if (m[i][j].score > highestScore) {
					highestScore = m[i][j].score;
					highestStruct = m[i][j];
				}
			}
		}
				
		vector<score_struct> path;
		
		int rowIndex = highestStruct.row;
		int pos = highestStruct.col;
		int score = highestStruct.score;
		
		while (pos >= 0 && score > 0) {
			score_struct temp = m[rowIndex][pos];
			score = temp.score;
			
			if (score > 0) {	path.push_back(temp);	}
			
			rowIndex = temp.prev;
			pos--;
		}

		reverse(path.begin(), path.end());
	
		return path;
		
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "extractHighestPath");
		exit(1);
	}
}
//***************************************************************************************************************
vector<trace_struct> Maligner::mapTraceRegionsToAlignment(vector<score_struct> path, vector<Sequence*> seqs) {
	try {
		vector<trace_struct> trace;
		
		int region_index = path[0].row;
		int region_start = path[0].col;
	
		for (int i = 1; i < path.size(); i++) {
		
			int next_region_index = path[i].row;
			
			if (next_region_index != region_index) {
				
				// add trace region
				int col_index = path[i].col;
				trace_struct temp;
				temp.col = region_start;
				temp.oldCol = col_index-1;
				temp.row = region_index;
				
				trace.push_back(temp);
							
				region_index = path[i].row;
				region_start = col_index;
			}
		}
	
		// get last one
		trace_struct temp;
		temp.col = region_start;
		temp.oldCol = path[path.size()-1].col;
		temp.row = region_index;
		trace.push_back(temp);

		return trace;
		
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "mapTraceRegionsToAlignment");
		exit(1);
	}
}
//***************************************************************************************************************
string Maligner::constructChimericSeq(vector<trace_struct> trace, vector<Sequence*> seqs) {
	try {
		string chimera = "";
		
		for (int i = 0; i < trace.size(); i++) {
			string seqAlign = seqs[trace[i].row]->getAligned();
			seqAlign = seqAlign.substr(trace[i].col, (trace[i].oldCol-trace[i].col+1));
			chimera += seqAlign;
		}
			
		return chimera;
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "constructChimericSeq");
		exit(1);
	}
}
//***************************************************************************************************************
float Maligner::computePercentID(string queryAlign, string chimera) {
	try {
	
		if (queryAlign.length() != chimera.length()) {
			mothurOut("Error, alignment strings are of different lengths: "); mothurOutEndLine();
			mothurOut(toString(queryAlign.length())); mothurOutEndLine(); mothurOutEndLine(); mothurOutEndLine(); mothurOutEndLine();
			mothurOut(toString(chimera.length())); mothurOutEndLine();
			return -1.0;
		}

	
		int numBases = 0;
		int numIdentical = 0;
	
		for (int i = 0; i < queryAlign.length(); i++) {
			if ((isalpha(queryAlign[i])) || (isalpha(chimera[i])))  {
				numBases++;		
				if (queryAlign[i] == chimera[i]) {
					numIdentical++;
				}
			}
		}
	
		if (numBases == 0) { return 0; }
	
		float percentIdentical = (numIdentical/(float)numBases) * 100;

		return percentIdentical;
		
	}
	catch(exception& e) {
		errorOut(e, "Maligner", "computePercentID");
		exit(1);
	}
}
//***************************************************************************************************************

