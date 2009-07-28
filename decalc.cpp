/*
 *  decalc.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "decalc.h"

//***************************************************************************************************************
void DeCalculator::setMask(string m) { 
	try {
		seqMask = m; 
		
		//whereever there is a base in the mask, save that value is query and subject
		for (int i = 0; i < seqMask.length(); i++) {
			if (isalpha(seqMask[i])) {
				h.insert(i);
			}
		}

	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "setMask");
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
		errorOut(e, "DeCalculator", "runMask");
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
		for (int i = 0; i < q.length(); i++) {
			if (isalpha(q[i]) && isalpha(s[i])) { front = i; break;  }
		}
		
		int back = 0;		
		for (int i = q.length(); i >= 0; i--) {
			if (isalpha(q[i]) && isalpha(s[i])) { back = i; break;  }
		}
		
		trim[front] = back;
		
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "trimSeqs");
		exit(1);
	}
}
//***************************************************************************************************************
//find the window breaks for each sequence - this is so you can move ahead by bases.
vector<int>  DeCalculator::findWindows(Sequence* query, int front, int back, int& size, int increment) {
	try {
		
		vector<int> win; 
		
		int cutoff = back - front;  //back - front 
			
		//if window is set to default
		if (size == 0) {  if (cutoff > 1200) {  size = 300; }
							else{  size = (cutoff / 4); }  } 
		else if (size > (cutoff / 4)) { 
				mothurOut("You have selected to large a window size for sequence " + query->getName() + ".  I will choose an appropriate window size."); mothurOutEndLine();
				size = (cutoff / 4); 
		}
	
		string seq = query->getAligned().substr(front, cutoff);
			
		//count bases
		int numBases = 0;
		for (int l = 0; l < seq.length(); l++) {  if (isalpha(seq[l])) { numBases++; }  }
			
		//save start of seq
		win.push_back(front);
		
		//move ahead increment bases at a time until all bases are in a window
		int countBases = 0;
		int totalBases = 0;  //used to eliminate window of blanks at end of sequence
			
		seq = query->getAligned();
		for (int m = front; m < (back - size) ; m++) {
				
			//count number of bases you see
			if (isalpha(seq[m])) { countBases++; totalBases++;  }
				
			//if you have seen enough bases to make a new window
			if (countBases >= increment) {
				win.push_back(m);  //save spot in alignment
				countBases = 0;				//reset bases you've seen in this window
			}
				
			//no need to continue if all your bases are in a window
			if (totalBases == numBases) {   break;  }
		}
			
		return win;
	
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "findWindows");
		exit(1);
	}
}

//***************************************************************************************************************
vector<float> DeCalculator::calcObserved(Sequence* query, Sequence* subject, vector<int> window, int size) {
	try {
		
		vector<float> temp;
//cout << "query length = " << query->getAligned().length() << '\t' << " subject length = " << subject.getAligned().length() << endl;				
		for (int m = 0; m < window.size(); m++) {
						
			string seqFrag = query->getAligned().substr(window[m], size);
			string seqFragsub = subject->getAligned().substr(window[m], size);
	//cout << "start point = " << window[m] << " end point = " << window[m]+size << endl;						
			int diff = 0;
			for (int b = 0; b < seqFrag.length(); b++) {
		
				if (seqFrag[b] != seqFragsub[b]) { diff++; }
			}
               
			//percentage of mismatched bases
			float dist;
			dist = diff / (float) seqFrag.length() * 100;       
				
			temp.push_back(dist);
		}
			
		return temp;
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "calcObserved");
		exit(1);
	}
}
//***************************************************************************************************************
float DeCalculator::calcDist(Sequence* query, Sequence* subject, int front, int back) {
	try {
		
		//so you only look at the trimmed part of the sequence
		int cutoff = back - front;  
			
		//from first startpoint with length back-front
		string seqFrag = query->getAligned().substr(front, cutoff);
		string seqFragsub = subject->getAligned().substr(front, cutoff);
														
		int diff = 0;
		for (int b = 0; b < seqFrag.length(); b++) {
			if (seqFrag[b] != seqFragsub[b]) { diff++; }
		}
               
		//percentage of mismatched bases
		float dist = diff / (float) seqFrag.length() * 100;       
				
		return dist;
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "calcDist");
		exit(1);
	}
}

//***************************************************************************************************************
vector<float> DeCalculator::calcExpected(vector<float> qav, float coef) {
	try {
		
		//for each window
		vector<float> queryExpected;
			
		for (int m = 0; m < qav.size(); m++) {		
				
			float expected = qav[m] * coef;
				
			queryExpected.push_back(expected);	
		}
			
		return queryExpected;
				
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "calcExpected");
		exit(1);
	}
}
//***************************************************************************************************************
float DeCalculator::calcDE(vector<float> obs, vector<float> exp) {
	try {
		
		//for each window
		float sum = 0.0;  //sum = sum from 1 to m of (oi-ei)^2
		for (int m = 0; m < obs.size(); m++) { 		sum += ((obs[m] - exp[m]) * (obs[m] - exp[m]));		}
			
		float de = sqrt((sum / (obs.size() - 1)));
			
		return de;
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "calcDE");
		exit(1);
	}
}

//***************************************************************************************************************

vector<float> DeCalculator::calcFreq(vector<Sequence*> seqs, string filename) {
	try {

		vector<float> prob;
		string freqfile = getRootName(filename) + "freq";
		ofstream outFreq;
		
		openOutputFile(freqfile, outFreq);
		
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
			for (int m = 0; m < freq.size(); m++) {   if (freq[m] > highest) {  highest = freq[m];  }		}
			
			float highFreq;
			//subtract gaps to "ignore them"
			if ( (seqs.size() - gaps) == 0 ) {  highFreq = 1.0;  }			
			else { highFreq = highest / (float) (seqs.size() - gaps);	 }
						
			float Pi;
			Pi =  (highFreq - 0.25) / 0.75; 
			
			//cannot have probability less than 0.
			if (Pi < 0) { Pi = 0.0; }
			
			//saves this for later
			outFreq << i+1 << '\t' << highFreq << endl;
			
			if (h.count(i) > 0) {
	cout << i+1 << '\t' << highFreq << endl;
				prob.push_back(Pi); 
			}
		}
		
		outFreq.close();
		
		return prob;
				
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "calcFreq");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float>  DeCalculator::findQav(vector<int> window, int size, vector<float> probabilityProfile) {
	try {
		vector<float>  averages; 
				
		//for each window find average
		for (int m = 0; m < window.size(); m++) {
				
			float average = 0.0;
				
			//while you are in the window for this sequence
			int count = 0;
			for (int j = window[m]; j < (window[m]+size); j++) {   
				
				//is this a spot that is included in the mask
				if (h.count(j) > 0) {
					average += probabilityProfile[j];
					count++;
				}
			}
				
			average = average / count;
	
			//save this windows average
			averages.push_back(average);
		}
				
		return averages;
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "findQav");
		exit(1);
	}
}

//***************************************************************************************************************
vector< vector<float> > DeCalculator::getQuantiles(vector<Sequence*> seqs, vector<int> windowSizesTemplate, int window, vector<float> probProfile, int increment, int start, int end) {
	try {
		vector< vector<float> > quan; 
		
		//percentage of mismatched pairs 1 to 100
		quan.resize(100);
		
		
		//for each sequence
		for(int i = start; i < end; i++){
		
			mothurOut("Processing template sequence " + toString(i)); mothurOutEndLine();
			Sequence* query = seqs[i];
			
			//compare to every other sequence in template
			for(int j = 0; j < i; j++){
				
				Sequence* subject = seqs[j];
				
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
				
				//dist-1 because vector indexes start at 0.
				quan[dist-1].push_back(de);
				
			}
		}

		return quan;
						
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "findQav");
		exit(1);
	}
}

//***************************************************************************************************************
float DeCalculator::getCoef(vector<float> obs, vector<float> qav) {
	try {
	
		//find average prob for this seqs windows
		float probAverage = 0.0;
		for (int j = 0; j < qav.size(); j++) {   probAverage += qav[j];	}
		probAverage = probAverage / (float) qav.size();
		
		//find observed average 
		float obsAverage = 0.0;
		for (int j = 0; j < obs.size(); j++) {   obsAverage += obs[j];	}
		obsAverage = obsAverage / (float) obs.size();
//cout << "sum ai / m = " << probAverage << endl;		
//cout << "sum oi / m = " << obsAverage << endl;
		float coef = obsAverage / probAverage;
						
		return coef;
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "getCoef");
		exit(1);
	}
}
//***************************************************************************************************************


