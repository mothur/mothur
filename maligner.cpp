/*
 *  maligner.cpp
 *  Mothur
 *
 *  Created by westcott on 9/23/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "maligner.h"

/***********************************************************************/ //int num, int match, int misMatch, , string mode, Database* dataLeft, Database* dataRight
Maligner::Maligner(vector<Sequence*> temp, int match, int misMatch, float div, int ms, int minCov) : db(temp), matchScore(match), misMatchPenalty(misMatch), minDivR(div), minSimilarity(ms), minCoverage(minCov) { 
			//numWanted(num),  , searchMethod(mode), databaseLeft(dataLeft), databaseRight(dataRight)
			
			m = MothurOut::getInstance(); 
			
}
/***********************************************************************/
string Maligner::getResults(Sequence* q, DeCalculator* decalc) {
	try {
		
		outputResults.clear();
		
		//make copy so trimming doesn't destroy query from calling class - remember to deallocate
		query = new Sequence(q->getName(), q->getAligned());
		
		string chimera;
		
		//copy refSeqs so that filter does not effect original
		for(int i = 0; i < db.size(); i++) {  
			Sequence* newSeq = new Sequence(db[i]->getName(), db[i]->getAligned());
			refSeqs.push_back(newSeq);
		}
		
		refSeqs = minCoverageFilter(refSeqs);
		
		if (refSeqs.size() < 2)  { 
			for (int i = 0; i < refSeqs.size(); i++) {  delete refSeqs[i];	}
			percentIdenticalQueryChimera = 0.0;
			return "unknown"; 
		}
		
		int chimeraPenalty = computeChimeraPenalty();
		
		//fills outputResults
		chimera = chimeraMaligner(chimeraPenalty, decalc);
		
		if (m->control_pressed) { return chimera;  }
				
		//free memory
		delete query;

		for (int i = 0; i < refSeqs.size(); i++) {  delete refSeqs[i];	}
		
		return chimera;
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "getResults");
		exit(1);
	}
}
/***********************************************************************/
string Maligner::chimeraMaligner(int chimeraPenalty, DeCalculator* decalc) {
	try {
		
		string chimera;
		
		//trims seqs to first non gap char in all seqs and last non gap char in all seqs
		spotMap = decalc->trimSeqs(query, refSeqs);
		
		//you trimmed the whole sequence, skip
		if (query->getAligned() == "") { return "no"; }

		vector<Sequence*> temp = refSeqs;
		temp.push_back(query);
			
		verticalFilter(temp);
		
		//for (int i = 0; i < refSeqs.size(); i++) { cout << refSeqs[i]->getName() << endl << refSeqs[i]->getAligned() << endl; }

		vector< vector<score_struct> > matrix = buildScoreMatrix(query->getAligned().length(), refSeqs.size()); //builds and initializes
		
		if (m->control_pressed) { return chimera;  }
		
		fillScoreMatrix(matrix, refSeqs, chimeraPenalty);
		
		vector<score_struct> path = extractHighestPath(matrix);
		
		if (m->control_pressed) { return chimera;  }
		
		vector<trace_struct> trace = mapTraceRegionsToAlignment(path, refSeqs);
		
		if (trace.size() > 1) {		chimera = "yes";	}
		else { chimera = "no";	return chimera; }
		
		int traceStart = path[0].col;
		int traceEnd = path[path.size()-1].col;	
		string queryInRange = query->getAligned();
		queryInRange = queryInRange.substr(traceStart, (traceEnd-traceStart));
		
		string chimeraSeq = constructChimericSeq(trace, refSeqs);
		
		percentIdenticalQueryChimera = computePercentID(queryInRange, chimeraSeq);
		
	/*	
		vector<trace_struct> trace = extractHighestPath(matrix);
				
		//cout << "traces\n";
		//for(int i=0;i<trace.size();i++){
		//	cout << trace[i].col << '\t' << trace[i].oldCol << '\t' << refSeqs[trace[i].row]->getName() << endl;
		//}
		
		if (trace.size() > 1) {		chimera = "yes";	}
		else { chimera = "no";	return chimera; }
		
		int traceStart = trace[0].col;
		int traceEnd = trace[trace.size()-1].oldCol;	
		string queryInRange = query->getAligned();
		queryInRange = queryInRange.substr(traceStart, (traceEnd-traceStart));*/
		
		if (m->control_pressed) { return chimera;  }
		
		//save output results
		for (int i = 0; i < trace.size(); i++) {
			int regionStart = trace[i].col;
			int regionEnd = trace[i].oldCol;
			int seqIndex = trace[i].row;
			
			results temp;
			
			temp.parent = refSeqs[seqIndex]->getName();
			temp.parentAligned = db[seqIndex]->getAligned();
			temp.nastRegionStart = spotMap[regionStart];
			temp.nastRegionEnd = spotMap[regionEnd];
			temp.regionStart = regionStart;
			temp.regionEnd = regionEnd;
			
			string parentInRange = refSeqs[seqIndex]->getAligned();
			parentInRange = parentInRange.substr(traceStart, (traceEnd-traceStart));
			
			temp.queryToParent = computePercentID(queryInRange, parentInRange);
			temp.divR = (percentIdenticalQueryChimera / temp.queryToParent);

			string queryInRegion = query->getAligned();
			queryInRegion = queryInRegion.substr(regionStart, (regionEnd-regionStart));
			
			string parentInRegion = refSeqs[seqIndex]->getAligned();
			parentInRegion = parentInRegion.substr(regionStart, (regionEnd-regionStart));
			
			temp.queryToParentLocal = computePercentID(queryInRegion, parentInRegion);
			
//			cout << temp.parent << '\t' << "NAST:" << temp.nastRegionStart << '-' << temp.nastRegionEnd << " G:" << temp.queryToParent << " L:" << temp.queryToParentLocal << endl;
			outputResults.push_back(temp);
		}
		
		return chimera;
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "chimeraMaligner");
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
			}else {
				delete ref[i];
			}
		}
		
		return newRefs;
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "minCoverageFilter");
		exit(1);
	}
}
/***********************************************************************/
// a breakpoint should yield fewer mismatches than this number with respect to the best parent sequence.
int Maligner::computeChimeraPenalty() {
	try {
		
		int numAllowable = ((1.0 - (1.0/minDivR)) * query->getNumBases());

//		if(numAllowable < 1){	numAllowable = 1;	}
		
		int penalty = int(numAllowable + 1) * misMatchPenalty;

		return penalty;

	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "computeChimeraPenalty");
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
		
		map<int, int> newMap;
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			string newAligned = "";
			int count = 0;
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is not a gap
				if (filterString[j] == '1') { 
					newAligned += seqAligned[j]; 
					newMap[count] = spotMap[j];
					count++;
				}
			}
			
			seqs[i]->setAligned(newAligned);
		}

		spotMap = newMap;
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "verticalFilter");
		exit(1);
	}
}
//***************************************************************************************************************
vector< vector<score_struct> > Maligner::buildScoreMatrix(int cols, int rows) {
	try{
		
		vector< vector<score_struct> > m(rows);
		
		for (int i = 0; i < rows; i++) {
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
		m->errorOut(e, "Maligner", "buildScoreMatrix");
		exit(1);
	}
}

//***************************************************************************************************************

void Maligner::fillScoreMatrix(vector<vector<score_struct> >& ms, vector<Sequence*> seqs, int penalty) {
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
				ms[i][0].score = 0;
//				ms[i][0].mismatches = 0;
			}else if (queryAligned[0] == subjectAligned[0])  { //|| subjectAligned[0] == 'N')
				ms[i][0].score = matchScore;
//				ms[i][0].mismatches = 0;
			}else{
				ms[i][0].score = 0;
//				ms[i][0].mismatches = 1;
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
					//matchMisMatchScore = matchScore;
					//leave the same
				}else if (queryAligned[j] == subjectAligned[j]) {
					matchMisMatchScore = matchScore;
//					ms[i][j].mismatches = ms[i][j-1].mismatches;
				}else if (queryAligned[j] != subjectAligned[j]) {
					matchMisMatchScore = misMatchPenalty;
//					ms[i][j].mismatches = ms[i][j-1].mismatches + 1;
				}
				
				//compute score based on previous columns scores
				for (int prevIndex = 0; prevIndex < numRows; prevIndex++) { //iterate through rows
					
					int sumScore = matchMisMatchScore + ms[prevIndex][j-1].score;
					
					//you are not at yourself
					if (prevIndex != i) {   sumScore += penalty;	}
					if (sumScore < 0)	{	sumScore = 0;			}
					
					if (sumScore > ms[i][j].score) {
						ms[i][j].score = sumScore;
						ms[i][j].prev = prevIndex;
					}
				}
			}
		}
		
	/*	for(int i=0;i<numRows;i++){
			cout << seqs[i]->getName();
			for(int j=0;j<numCols;j++){
				cout << '\t' << ms[i][j].mismatches;
			}
			cout << endl;
		}
		cout << endl;*/
		/*cout << numRows << '\t' << numCols << endl;
		for(int i=0;i<numRows;i++){
			cout << seqs[i]->getName() << endl << seqs[i]->getAligned() << endl << endl;
			if ((seqs[i]->getName() == "S000003470") || (seqs[i]->getName() == "S000383265") || (seqs[i]->getName() == "7000004128191054")) {
			for(int j=0;j<numCols;j++){
				cout << '\t' << ms[i][j].score;
			}
			cout << endl;
			}
		}
		cout << endl;*/
		/*for(int i=0;i<numRows;i++){
			cout << seqs[i]->getName();
			for(int j=0;j<numCols;j++){
				cout << '\t' << ms[i][j].prev;
			}
			cout << endl;
		}*/
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "fillScoreMatrix");
		exit(1);
	}
}
//***************************************************************************************************************
vector<score_struct> Maligner::extractHighestPath(vector<vector<score_struct> > ms) {
	try {
		
		//get matrix dimensions
		int numCols = query->getAligned().length();
		int numRows = ms.size();
		
		
		//find highest score scoring matrix
		score_struct highestStruct;
		int highestScore = 0;
		
		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				if (ms[i][j].score > highestScore) {
					highestScore = ms[i][j].score;
					highestStruct = ms[i][j];
				}
			}
		}
		
		vector<score_struct> path;
		
		int rowIndex = highestStruct.row;
		int pos = highestStruct.col;
		int score = highestStruct.score;
		
		while (pos >= 0 && score > 0) {
			score_struct temp = ms[rowIndex][pos];
			score = temp.score;
			
			if (score > 0) {	path.push_back(temp);	}
			
			rowIndex = temp.prev;
			pos--;
		}
		
		reverse(path.begin(), path.end());
		
		return path;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "extractHighestPath");
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
		m->errorOut(e, "Maligner", "mapTraceRegionsToAlignment");
		exit(1);
	}
}

/***************************************************************************************************************

vector<trace_struct> Maligner::extractHighestPath(vector<vector<score_struct> > ms) {
	try {
	
		
		//get matrix dimensions
		int numCols = query->getAligned().length();
		int numRows = ms.size();
	
	
		//find highest score scoring matrix
		vector<score_struct> highestStruct;
		int highestScore = 0;
		
		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				if (ms[i][j].score > highestScore) {
					highestScore = ms[i][j].score;
					highestStruct.resize(0);
					highestStruct.push_back(ms[i][j]);
				}
				else if(ms[i][j].score == highestScore){
					highestStruct.push_back(ms[i][j]);
				}
			}
		}
			
		//cout << endl << highestScore << '\t' << highestStruct.size() << '\t' << highestStruct[0].row << endl;	
		
		vector<trace_struct> maxTrace;
		double maxPercentIdenticalQueryAntiChimera = 0;
		
		for(int i=0;i<highestStruct.size();i++){
			
			vector<score_struct> path;

			int rowIndex = highestStruct[i].row;
			int pos = highestStruct[i].col;
			int score = highestStruct[i].score;
					
			while (pos >= 0 && score > 0) {
				score_struct temp = ms[rowIndex][pos];
				score = temp.score;
				
				if (score > 0) {	path.push_back(temp);	}
				
				rowIndex = temp.prev;
				pos--;
			}

			reverse(path.begin(), path.end());

			vector<trace_struct> trace = mapTraceRegionsToAlignment(path, refSeqs);
		
			//cout << "traces\n";
			//for(int j=0;j<trace.size();j++){
			//	cout << trace[j].col << '\t' << trace[j].oldCol << '\t' << refSeqs[trace[j].row]->getName() << endl;
			//}
						
			int traceStart = path[0].col;
			int traceEnd = path[path.size()-1].col;	
//			cout << "traceStart/End\t" << traceStart << '\t' << traceEnd << endl;
			
			string queryInRange = query->getAligned();
			queryInRange = queryInRange.substr(traceStart, (traceEnd-traceStart+1));
//			cout << "here" << endl;
			string chimeraSeq = constructChimericSeq(trace, refSeqs);
			string antiChimeraSeq = constructAntiChimericSeq(trace, refSeqs);
		
			percentIdenticalQueryChimera = computePercentID(queryInRange, chimeraSeq);			
			double percentIdenticalQueryAntiChimera = computePercentID(queryInRange, antiChimeraSeq);
//			cout << i << '\t' << percentIdenticalQueryChimera << '\t' << percentIdenticalQueryAntiChimera << endl;
			
			if(percentIdenticalQueryAntiChimera > maxPercentIdenticalQueryAntiChimera){
				maxPercentIdenticalQueryAntiChimera = percentIdenticalQueryAntiChimera;
				maxTrace = trace;
			}
		}
//		cout << maxPercentIdenticalQueryAntiChimera << endl;
		return maxTrace;
	
		
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "extractHighestPath");
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

//		cout << endl;
//		cout << trace.size() << endl;
//		for(int i=0;i<trace.size();i++){
//			cout << seqs[trace[i].row]->getName() << endl;
//		}
//		cout << endl;
		
		return trace;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "mapTraceRegionsToAlignment");
		exit(1);
	}
}
*/
//***************************************************************************************************************

string Maligner::constructChimericSeq(vector<trace_struct> trace, vector<Sequence*> seqs) {
	try {
		string chimera = "";
		
		for (int i = 0; i < trace.size(); i++) {
//			cout << i << '\t' << trace[i].row << '\t' << trace[i].col << '\t' << trace[i].oldCol << endl;
			
			string seqAlign = seqs[trace[i].row]->getAligned();
			seqAlign = seqAlign.substr(trace[i].col, (trace[i].oldCol-trace[i].col+1));
			chimera += seqAlign;
		}
		
		if (chimera != "") { chimera = chimera.substr(0, (chimera.length()-1)); }
		return chimera;
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "constructChimericSeq");
		exit(1);
	}
}

//***************************************************************************************************************

string Maligner::constructAntiChimericSeq(vector<trace_struct> trace, vector<Sequence*> seqs) {
	try {
		string antiChimera = "";
		
		for (int i = 0; i < trace.size(); i++) {
//			cout << i << '\t' << (trace.size() - i - 1) << '\t' << trace[i].row << '\t' << trace[i].col << '\t' << trace[i].oldCol << endl;
			
			int oppositeIndex = trace.size() - i - 1;
			
			string seqAlign = seqs[trace[oppositeIndex].row]->getAligned();
			seqAlign = seqAlign.substr(trace[i].col, (trace[i].oldCol-trace[i].col+1));
			antiChimera += seqAlign;
		}
		
		return antiChimera;
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "constructChimericSeq");
		exit(1);
	}
}

//***************************************************************************************************************
float Maligner::computePercentID(string queryAlign, string chimera) {
	try {
	
		if (queryAlign.length() != chimera.length()) {
			m->mothurOut("Error, alignment strings are of different lengths: "); m->mothurOutEndLine();
			m->mothurOut(toString(queryAlign.length())); m->mothurOutEndLine(); 
			m->mothurOut(toString(chimera.length())); m->mothurOutEndLine();
			return -1.0;
		}

		int numIdentical = 0;
		int countA = 0;
		int countB = 0;
		for (int i = 0; i < queryAlign.length(); i++) {
			if (((queryAlign[i] != 'G') && (queryAlign[i] != 'T') && (queryAlign[i] != 'A') && (queryAlign[i] != 'C')&& (queryAlign[i] != '.') && (queryAlign[i] != '-')) ||
				((chimera[i] != 'G') && (chimera[i] != 'T') && (chimera[i] != 'A') && (chimera[i] != 'C')&& (chimera[i] != '.') && (chimera[i] != '-'))) {}
			else {

				bool charA = false; bool charB = false;
				if ((queryAlign[i] == 'G') || (queryAlign[i] == 'T') || (queryAlign[i] == 'A') || (queryAlign[i] == 'C')) { charA = true; }
				if ((chimera[i] == 'G') || (chimera[i] == 'T') || (chimera[i] == 'A') || (chimera[i] == 'C')) { charB = true; }
					
				if (charA || charB) {
						
					if (charA) { countA++; }
					if (charB) { countB++; }
						
					if (queryAlign[i] == chimera[i]) {
						numIdentical++;
					}
				}
			}
		}
		
		float numBases = (countA + countB) /(float) 2;
		
		if (numBases == 0) { return 0; }
	
		float percentIdentical = (numIdentical/(float)numBases) * 100;

		return percentIdentical;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Maligner", "computePercentID");
		exit(1);
	}
}
//***************************************************************************************************************
