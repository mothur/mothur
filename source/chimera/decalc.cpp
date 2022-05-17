/*
 *  decalc.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "decalc.h"
#include "mothurchimera.h"
#include "calculator.h"
#include "eachgapdist.h"
#include "ignoregaps.h"
#include "eachgapdist.h"

//***************************************************************************************************************
void DeCalculator::setMask(string ms) { 
	try {
		seqMask = ms; 
		int count = 0;
		maskMap.clear();
		
		if (seqMask.length() != 0) {
			//whereever there is a base in the mask, save that value is query and subject
			for (int i = 0; i < seqMask.length(); i++) {
				if (isalpha(seqMask[i])) {
					h.insert(i);
					maskMap[count] = i;
					count++;
					
				}
			}
		}else {
			for (int i = 0; i < alignLength; i++) {   
				h.insert(i); 
				maskMap[count] = i;
				count++;
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "setMask");
		exit(1);
	} 
}
//***************************************************************************************************************
void DeCalculator::runMask(Sequence* seq) {
	try{
		
		string q = seq->getAligned();
		string tempQuery = "";
		
		//whereever there is a base in the mask, save that value is query and subject
		set<int>::iterator setit;
		for ( setit=h.begin() ; setit != h.end(); setit++ )  {
			tempQuery += q[*setit];
		}
		
		//save masked values
		seq->setAligned(tempQuery);
		seq->setUnaligned(tempQuery);
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "runMask");
		exit(1);
	}
}
//***************************************************************************************************************
//num is query's spot in querySeqs
void DeCalculator::trimSeqs(Sequence* query, Sequence* subject, map<int, int>& trim) {
	try {
		string q = query->getAligned();
		string s = subject->getAligned();
		
		int front = 0;
		for (int i = 0; i < q.length(); i++) { if (isalpha(q[i]) && isalpha(s[i])) { front = i; break;  } }
	
		int back = 0;		
		for (int i = q.length(); i >= 0; i--) { if (isalpha(q[i]) && isalpha(s[i])) { back = i; break;  } }
		
		trim[front] = back;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "trimSeqs");
		exit(1);
	}
}
//***************************************************************************************************************
vector<int>  DeCalculator::findWindows(Sequence* query, int front, int back, int& size, int increment) {
	try {
		vector<int> win; 
		
		int cutoff = back - front;  //back - front 
			
		//if window is set to default
		if (size == 0) {  if (cutoff > 1200) {  size = 300; }
							else{  size = (cutoff / 4); }  } 
		else if (size > (cutoff / 4)) { 
				m->mothurOut("[WARNING]: You have selected too large a window size for sequence " + query->getName() + ".  I will choose an appropriate window size.\n");
				size = (cutoff / 4); 
		}
	
		//this follows wigeon, but we may want to consider that it chops off the end values if the sequence cannot be evenly divided into steps
		for (int i = front; i < (back - size) ; i+=increment) {  win.push_back(i);  }
        
		return win;
	
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "findWindows");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float> DeCalculator::calcObserved(Sequence* query, Sequence* subject, vector<int> window, int size) {
	try {
		
        if (query->getAligned().length() != subject->getAligned().length()) {
            
            m->mothurOut("[ERROR]: seqLengths must match. queryLength = " + toString(query->getAligned().length()) + ", subjectLength = " + toString(subject->getAligned().length()) +"\n");
            m->setControl_pressed(true);
        }
        
        
		vector<float> temp;	
		//int gaps = 0;		
		for (int i = 0; i < window.size(); i++) {
					
            if (m->getControl_pressed()) { break; }
            
			string seqFrag = query->getAligned().substr(window[i], size);
			string seqFragsub = subject->getAligned().substr(window[i], size);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: seqLengths = " + toString(seqFrag.length()) + "," + toString(seqFragsub.length()) +  "\t" + toString(window[i]) + "\t" + toString(size) +"\n"); }
				
			int diff = 0;
			for (int b = 0; b < seqFrag.length(); b++) {
				//if at least one is a base and they are not equal
				if( (isalpha(seqFrag[b]) || isalpha(seqFragsub[b])) && (seqFrag[b] != seqFragsub[b]) ) { diff++; }
            }
            
			//percentage of mismatched bases
			float dist = diff / (float) (seqFrag.length()) * 100;
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: diffs = " + toString(diff) + ", dist = " + toString(dist) + "\n"); }
			
			temp.push_back(dist);
		}
			
		return temp;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "calcObserved");
		exit(1);
	}
}
//***************************************************************************************************************
float DeCalculator::calcDist(Sequence* query, Sequence* subject, int front, int back) {
	try {
		
		//so you only look at the trimmed part of the sequence
		int cutoff = back - front;  
		int gaps = 0;
			
		//from first startpoint with length back-front
		string seqFrag = query->getAligned().substr(front, cutoff);
		string seqFragsub = subject->getAligned().substr(front, cutoff);
														
		int diff = 0;
		for (int b = 0; b < seqFrag.length(); b++) {
			//ignore gaps
			if((!isalpha(seqFrag[b])) && (!isalpha(seqFragsub[b]))) { gaps++; }
			if (seqFrag[b] != seqFragsub[b]) { diff++; }
		}
		
		//if the whole fragment is 0 distance = 0
		if ((seqFrag.length()-gaps) == 0)  { return 0.0; }
               
		//percentage of mismatched bases
		float dist = diff / (float) (seqFrag.length()-gaps) * 100;       
				
		return dist;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "calcDist");
		exit(1);
	}
}

//***************************************************************************************************************
vector<float> DeCalculator::calcExpected(vector<float> qav, float coef) {
	try {
		
		//for each window
		vector<float> queryExpected;
			
		for (int j = 0; j < qav.size(); j++) {		
				
			float expected = qav[j] * coef;
				
			queryExpected.push_back(expected);	
		}
			
		return queryExpected;
				
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "calcExpected");
		exit(1);
	}
}
//***************************************************************************************************************
float DeCalculator::calcDE(vector<float> obs, vector<float> exp) {
	try {
		
		//for each window
		float sum = 0.0;  //sum = sum from 1 to i of (oi-ei)^2
		int numZeros = 0;
		for (int j = 0; j < obs.size(); j++) { 	sum += ((obs[j] - exp[j]) * (obs[j] - exp[j]));	}
			
		float de = sqrt((sum / (obs.size() - 1 - numZeros)));
			
		return de;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "calcDE");
		exit(1);
	}
}

//***************************************************************************************************************

vector<float> DeCalculator::calcFreq(vector<Sequence*> seqs, string filename, string version) {
	try {

		vector<float> prob;
        Utils util;
		string freqfile = util.getRootName(filename) + "freq";
		ofstream outFreq;
		
		util.openOutputFile(freqfile, outFreq);
		
		outFreq << "#" << version << endl;
		
		string length = toString(seqs.size());  //if there are 5000 seqs in the template then set precision to 3
		int precision = length.length() - 1;
		
		//format output
		outFreq.setf(ios::fixed, ios::floatfield); outFreq.setf(ios::showpoint);
		
		//at each position in the sequence
		for (int i = 0; i < seqs[0]->getAligned().length(); i++) {
			
			vector<int> freq;   freq.resize(4,0);
			int gaps = 0;
			
			//find the frequency of each nucleotide
			for (int j = 0; j < seqs.size(); j++) {
				
				char value = seqs[j]->getAligned()[i];
				
				if(toupper(value) == 'A')									{	freq[0]++;	}
				else if(toupper(value) == 'T' || toupper(value) == 'U')		{	freq[1]++;	}
				else if(toupper(value) == 'G')								{	freq[2]++;	}
				else if(toupper(value) == 'C')								{	freq[3]++;	}
				else { gaps++; }
			}
			
			//find base with highest frequency
			int highest = 0;
			for (int j = 0; j < freq.size(); j++) {   if (freq[j] > highest) {  highest = freq[j];  }		}
			
			float highFreq = highest / (float) (seqs.size());	
			
			float Pi;
			Pi =  (highFreq - 0.25) / 0.75; 
			
			//cannot have probability less than 0.
			if (Pi < 0) { Pi = 0.0; }
			
			//saves this for later
			outFreq << setprecision(precision) << i << '\t' << highFreq << endl;
	
			if (h.count(i) > 0) {
				prob.push_back(Pi); 
			}
		}
		
		outFreq.close();
		
		return prob;
				
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "calcFreq");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float>  DeCalculator::findQav(vector<int> window, int size, vector<float> probabilityProfile) {
	try {
		vector<float>  averages; 
				
		//for each window find average
		for (int i = 0; i < window.size(); i++) {
				
			float average = 0.0;
				
			//while you are in the window for this sequence
			int count = 0;
			for (int j = window[i]; j < (window[i]+size); j++) {   
				average += probabilityProfile[j];
				count++;
			}
				
			average = average / count;
	
			//save this windows average
			averages.push_back(average);
		}
				
		return averages;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "findQav");
		exit(1);
	}
}
//***************************************************************************************************************
//seqs have already been masked
vector< vector<float> > DeCalculator::getQuantiles(vector<Sequence*> seqs, vector<int> windowSizesTemplate, int window, vector<float> probProfile, int increment, int start, int end) {
	try {
		vector< vector<float> > quan; 
		
		//percentage of mismatched pairs 1 to 100
		quan.resize(100);

		//for each sequence
		for(int i = start; i < end; i++){
		
			m->mothurOut("Processing sequence " + toString(i)+ "\n");
			Sequence* query = new Sequence(seqs[i]->getName(), seqs[i]->getAligned());
			
			//compare to every other sequence in template
			for(int j = 0; j < i; j++){
				
				Sequence* subject = new Sequence(seqs[j]->getName(), seqs[j]->getAligned());
				
				if (m->getControl_pressed()) { delete query; delete subject; return quan; }
				
				map<int, int> trim;
				map<int, int>::iterator it;
				
				trimSeqs(query, subject, trim);
				
				it = trim.begin();
				int front = it->first; int back = it->second;
				
				//reset window for each new comparison
				windowSizesTemplate[i] = window;
				
				vector<int> win = findWindows(query, front, back, windowSizesTemplate[i], increment);
				
				vector<float> obsi = calcObserved(query, subject, win, windowSizesTemplate[i]);
				
				vector<float> q = findQav(win, windowSizesTemplate[i], probProfile);
									
				float alpha = getCoef(obsi, q);
						
				vector<float> exp = calcExpected(q, alpha);
				
				float de = calcDE(obsi, exp);
								
				float dist = calcDist(query, subject, front, back); 
				
				dist = ceil(dist);

                quan[dist].push_back(de);

				delete subject;
			}
			
			delete query;
		}
	

		return quan;
						
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "getQuantiles");
		exit(1);
	}
}
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareQuanMembers(quanMember left, quanMember right){
	return (left.score < right.score);	
} 
//***************************************************************************************************************
//this was going to be used by pintail to increase the sensitivity of the chimera detection, but it wasn't quite right.  may want to revisit in the future...
void DeCalculator::removeObviousOutliers(vector< vector<float> >& quantiles, int num) {
	try {
						
		for (int i = 0; i < quantiles.size(); i++) {
			
			//find mean of this quantile score
			sort(quantiles[i].begin(), quantiles[i].end());
			
			vector<float> temp;
			if (quantiles[i].size() != 0) {
				float high = quantiles[i][int(quantiles[i].size() * 0.99)];
				float low =  quantiles[i][int(quantiles[i].size() * 0.01)];
			
				//look at each value in quantiles to see if it is an outlier
				for (int j = 0; j < quantiles[i].size(); j++) {
					//is this score between 1 and 99%
					if ((quantiles[i][j] > low) && (quantiles[i][j] < high)) {
						temp.push_back(quantiles[i][j]);
					}
				}
			}
			quantiles[i] = temp;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "removeObviousOutliers");
		exit(1);
	}
}
//***************************************************************************************************************
float DeCalculator::findAverage(vector<float> myVector) {
	try{
		
		float total = 0.0;
		for (int i = 0; i < myVector.size(); i++)  {  total += myVector[i];  }
		
		float average = total / (float) myVector.size();
		
		return average;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "findAverage");
		exit(1);
	}
}

//***************************************************************************************************************
float DeCalculator::getCoef(vector<float> obs, vector<float> qav) {
	try {
	
		//find average prob for this seqs windows
		float probAverage = findAverage(qav);
				
		//find observed average 
		float obsAverage = findAverage(obs);
			
		float coef = obsAverage / probAverage;
						
		return coef;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "getCoef");
		exit(1);
	}
}
//***************************************************************************************************************
//gets closest matches to each end, since chimeras will most likely have different parents on each end
vector<Sequence> DeCalculator::findClosest(Sequence querySeq, vector<Sequence*>& thisTemplate, vector<Sequence*>& thisFilteredTemplate, int numWanted, int minSim) {
	try {
		//indexes.clear();
		
		vector<Sequence> seqsMatches;  
		
		vector<SeqDist> distsLeft;
		vector<SeqDist> distsRight;
		
        DistCalc* distcalculator = new eachGapDist(1.0);
		
		string queryUnAligned = querySeq.getUnaligned();
		int numBases = int(queryUnAligned.length() * 0.33);
		
		string leftQuery = ""; //first 1/3 of the sequence
		string rightQuery = ""; //last 1/3 of the sequence
		string queryAligned = querySeq.getAligned();
		
		//left side
		bool foundFirstBase = false;
		int baseCount = 0;
		int leftSpot = 0;
		int firstBaseSpot = 0;
		for (int i = 0; i < queryAligned.length(); i++) {
			//if you are a base
			if (isalpha(queryAligned[i])) {		
				baseCount++; 
				if (!foundFirstBase) {   foundFirstBase = true;  firstBaseSpot = i;  }
			}
			
			//eliminate opening .'s
			if (foundFirstBase) {   leftQuery += queryAligned[i];  }
			//if you have 1/3
			if (baseCount >= numBases) {  leftSpot = i; break; } //first 1/3
		}
		
		//right side - count through another 1/3, so you are at last third
		baseCount = 0;
		int rightSpot = 0;
		for (int i = leftSpot; i < queryAligned.length(); i++) {
			//if you are a base
			if (isalpha(queryAligned[i])) {		baseCount++;	}
			//if you have 1/3
			if (baseCount > numBases + 1) { rightSpot = i;  break; } //last 1/3
		}
		
		//trim end
		//find last position in query that is a non gap character
		int lastBaseSpot = queryAligned.length()-1;
		for (int j = queryAligned.length()-1; j >= 0; j--) {
			if (isalpha(queryAligned[j])) {
				lastBaseSpot = j;
				break;
			}
		}
		rightQuery = queryAligned.substr(rightSpot, (lastBaseSpot-rightSpot+1)); //sequence from pos spot to end
		
		Sequence queryLeft(querySeq.getName(), leftQuery);
		Sequence queryRight(querySeq.getName(), rightQuery);
		
		for(int j = 0; j < thisFilteredTemplate.size(); j++){
			
			string dbAligned = thisFilteredTemplate[j]->getAligned();
			string leftDB = dbAligned.substr(firstBaseSpot, (leftSpot-firstBaseSpot+1)); //first 1/3 of the sequence
			string rightDB = dbAligned.substr(rightSpot, (lastBaseSpot-rightSpot+1)); //last 1/3 of the sequence

			Sequence dbLeft(thisFilteredTemplate[j]->getName(), leftDB);
			Sequence dbRight(thisFilteredTemplate[j]->getName(), rightDB);

			double distLeft = distcalculator->calcDist(queryLeft, dbLeft);
			double distRight = distcalculator->calcDist(queryRight, dbRight);

			SeqDist subjectLeft;
			subjectLeft.seq = nullptr;
			subjectLeft.dist = distLeft;
			subjectLeft.index = j;
			
			distsLeft.push_back(subjectLeft);
			
			SeqDist subjectRight;
			subjectRight.seq = nullptr;
			subjectRight.dist = distRight;
			subjectRight.index = j;
			
			distsRight.push_back(subjectRight);

		}
		
		delete distcalculator;
		
		//sort by smallest distance
		sort(distsRight.begin(), distsRight.end(), compareSeqDist);
		sort(distsLeft.begin(), distsLeft.end(), compareSeqDist);

		
		//merge results		
		map<string, string> seen;
		map<string, string>::iterator it;
		
		vector<SeqDist> dists;
		float lastRight = distsRight[0].dist;
		float lastLeft = distsLeft[0].dist;

		float maxDist = 1.0 - (minSim / 100.0);

		for (int i = 0; i < numWanted+1; i++) {
			if (m->getControl_pressed()) { return seqsMatches; }
			
			//add left if you havent already
			it = seen.find(thisTemplate[distsLeft[i].index]->getName());
			if (it == seen.end() && distsLeft[i].dist <= maxDist) {  
				dists.push_back(distsLeft[i]);
				seen[thisTemplate[distsLeft[i].index]->getName()] = thisTemplate[distsLeft[i].index]->getName();
				lastLeft =  distsLeft[i].dist;
			}

			//add right if you havent already
			it = seen.find(thisTemplate[distsRight[i].index]->getName());
			if (it == seen.end() && distsRight[i].dist <= maxDist) {  
				dists.push_back(distsRight[i]);
				seen[thisTemplate[distsRight[i].index]->getName()] = thisTemplate[distsRight[i].index]->getName();
				lastRight =  distsRight[i].dist;
			}
			
			if (i == numWanted) { break; }
			
		}
		
		//are we still above the minimum similarity cutoff
		if ((lastLeft >= minSim) || (lastRight >= minSim)) {
			//add in ties from left
			int i = numWanted;
			while (i < distsLeft.size()) {  
				if (distsLeft[i].dist == lastLeft) {  dists.push_back(distsLeft[i]);  }
				else { break; }
				i++;
			}
			
			//add in ties from right
			i = numWanted;
			while (i < distsRight.size()) {  
				if (distsRight[i].dist == lastRight) {  dists.push_back(distsRight[i]);  }
				else { break; }
				i++;
			}
		}
		
		for (int i = 0; i < dists.size(); i++) {

			if ((thisTemplate[dists[i].index]->getName() != querySeq.getName()) && (((1.0-dists[i].dist)*100) >= minSim)) {
				Sequence temp(thisTemplate[dists[i].index]->getName(), thisTemplate[dists[i].index]->getAligned()); //have to make a copy so you can trim and filter without stepping on eachother.
				seqsMatches.push_back(temp);
			}

		}
		
		return seqsMatches;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "findClosest");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence* DeCalculator::findClosest(Sequence* querySeq, vector<Sequence*> db) {
	try {
		
		Sequence* seqsMatch;  
		
        DistCalc* distcalculator = new eachGapDist(1.0);
		int index = 0;
		int smallest = 1000000;
		
		for(int j = 0; j < db.size(); j++){
			
			double dist = distcalculator->calcDist(*querySeq, *db[j]);
			
			if (dist < smallest) { 
				smallest = dist;
				index = j;
			}
		}
		
		delete distcalculator;
		
		seqsMatch = new Sequence(db[index]->getName(), db[index]->getAligned()); //have to make a copy so you can trim and filter without stepping on eachother.
			
		return seqsMatch;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "findClosest");
		exit(1);
	}
}
/***************************************************************************************************************/
map<int, int> DeCalculator::trimSeqs(Sequence& query, vector<Sequence>& topMatches) {
	try {
		
		int frontPos = 0;  //should contain first position in all seqs that is not a gap character
		int rearPos = query.getAligned().length();
		
		//********find first position in topMatches that is a non gap character***********//
		//find first position all query seqs that is a non gap character
		for (int i = 0; i < topMatches.size(); i++) {
			
			string aligned = topMatches[i].getAligned();
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
		
		
		string aligned = query.getAligned();
		int pos = 0;
			
		//find first position in query that is a non gap character
		for (int j = 0; j < aligned.length(); j++) {
			if (isalpha(aligned[j])) {
				pos = j;
				break;
			}
		}
		
		//save this spot if it is the farthest
		if (pos > frontPos) { frontPos = pos; }
		
		
		//********find last position in topMatches that is a non gap character***********//
		for (int i = 0; i < topMatches.size(); i++) {
			
			string aligned = topMatches[i].getAligned();
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
		
		
		aligned = query.getAligned();
		pos = aligned.length();
		
		//find last position in query that is a non gap character
		for (int j = aligned.length()-1; j >= 0; j--) {
			if (isalpha(aligned[j])) {
				pos = j;
				break;
			}
		}
		
		//save this spot if it is the farthest
		if (pos < rearPos) { rearPos = pos; }
		
		map<int, int> trimmedPos;
		//check to make sure that is not whole seq
		if ((rearPos - frontPos - 1) <= 0) {  
			query.setAligned("");
			//trim topMatches
			for (int i = 0; i < topMatches.size(); i++) {
				topMatches[i].setAligned("");
			}
			
		}else {

			//trim query
			string newAligned = query.getAligned();
			newAligned = newAligned.substr(frontPos, (rearPos-frontPos+1));
			query.setAligned(newAligned);
			
			//trim topMatches
			for (int i = 0; i < topMatches.size(); i++) {
				newAligned = topMatches[i].getAligned();
				newAligned = newAligned.substr(frontPos, (rearPos-frontPos+1));
				topMatches[i].setAligned(newAligned);
			}
			
			for (int i = 0; i < newAligned.length(); i++) {
				trimmedPos[i] = i+frontPos;
			}
		}
		return trimmedPos;
	}
	catch(exception& e) {
		m->errorOut(e, "DeCalculator", "trimSequences");
		exit(1);
	}

}
//***************************************************************************************************************



