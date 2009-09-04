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
//cout << "query = " << q[i] << " subject = " << s[i] << endl;
			if (isalpha(q[i]) && isalpha(s[i])) { front = i; break;  }
		}
//cout << endl << endl;		
		int back = 0;		
		for (int i = q.length(); i >= 0; i--) {
//cout << "query = " << q[i] << " subject = " << s[i] << endl;
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
	
	/*	string seq = query->getAligned().substr(front, cutoff);
			
		//count bases
		int numBases = 0;
		for (int l = 0; l < seq.length(); l++) {  if (isalpha(seq[l])) { numBases++; }  }
//cout << "num Bases = " << numBases << endl;			
		//save start of seq
		win.push_back(front);
//cout << front << '\t';		
		//move ahead increment bases at a time until all bases are in a window
		int countBases = 0;
		int totalBases = 0;  //used to eliminate window of blanks at end of sequence
			
		seq = query->getAligned();
		for (int m = front; m < (back - size) ; m++) {
				
			//count number of bases you see
			if (isalpha(seq[m])) { countBases++;  }
				
			//if you have seen enough bases to make a new window
			if (countBases >= increment) {
				//total bases is the number of bases in a window already.
				totalBases += countBases;
//cout << "total bases = " << totalBases << endl;
				win.push_back(m);  //save spot in alignment
//cout << m << '\t';
				countBases = 0;				//reset bases you've seen in this window
			}
				
			//no need to continue if all your bases are in a window
			if (totalBases == numBases) {   break;  }
		}
	

		//get last window if needed
		if (totalBases < numBases) {   win.push_back(back-size);  }
//cout << endl << endl;
*/	
		
		//this follows wigeon, but we may want to consider that it chops off the end values if the sequence cannot be evenly divided into steps
		for (int m = front; m < (back - size) ; m+=increment) {  win.push_back(m);  }


	
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
		//int gaps = 0;		
		for (int m = 0; m < window.size(); m++) {
						
			string seqFrag = query->getAligned().substr(window[m], size);
			string seqFragsub = subject->getAligned().substr(window[m], size);
				
			int diff = 0;
			for (int b = 0; b < seqFrag.length(); b++) {
				//if at least one is a base and they are not equal
				if( (isalpha(seqFrag[b]) || isalpha(seqFragsub[b])) && (seqFrag[b] != seqFragsub[b]) ) { diff++; }
				
				//ignore gaps
				//if((!isalpha(seqFrag[b])) && (!isalpha(seqFragsub[b]))) { gaps++; }
			}
               
			//percentage of mismatched bases
			float dist;
			
			//if the whole fragment is 0 distance = 0
			//if ((seqFrag.length()-gaps) == 0)  { dist =  0.0; }
               
			//percentage of mismatched bases
			//else {  dist = diff / (float) (seqFrag.length()-gaps) * 100;   } 
			
			dist = diff / (float) (seqFrag.length()) * 100;	
			
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
		int numZeros = 0;
		for (int m = 0; m < obs.size(); m++) { 		
			
			//if (obs[m] != 0.0) {
				sum += ((obs[m] - exp[m]) * (obs[m] - exp[m]));	
			//}else {  numZeros++;   }
			
		}
			
		float de = sqrt((sum / (obs.size() - 1 - numZeros)));
			
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
			for (int m = 0; m < freq.size(); m++) {   if (freq[m] > highest) {  highest = freq[m];  }		}
			
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
		errorOut(e, "DeCalculator", "findQav");
		exit(1);
	}
}
//***************************************************************************************************************
//seqs have already been masked
vector< vector<quanMember> > DeCalculator::getQuantiles(vector<Sequence*> seqs, vector<int> windowSizesTemplate, int window, vector<float> probProfile, int increment, int start, int end) {
	try {
		vector< vector<quanMember> > quan; 
		
		//percentage of mismatched pairs 1 to 100
		quan.resize(100);
//ofstream o;
//string out = "getQuantiles.out";
//openOutputFile(out, o);
		
		//for each sequence
		for(int i = start; i < end; i++){
		
			mothurOut("Processing sequence " + toString(i)); mothurOutEndLine();
			Sequence* query = new Sequence(seqs[i]->getName(), seqs[i]->getAligned());
			
			//compare to every other sequence in template
			for(int j = 0; j < i; j++){
				
				Sequence* subject = new Sequence(seqs[j]->getName(), seqs[j]->getAligned());
				
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
	//o << i << '\t' <<  j << '\t' << dist << '\t' << de << endl;			
				dist = ceil(dist);
				
				quanMember newScore(de, i, j);
				
				//dist-1 because vector indexes start at 0.
				quan[dist-1].push_back(newScore);
				
				delete subject;
			}
			
			delete query;
		}

		return quan;
						
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "getQuantiles");
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
void DeCalculator::removeObviousOutliers(vector< vector<quanMember> >& quantiles, int num) {
	try {
						
		for (int i = 0; i < quantiles.size(); i++) {
		
			//find mean of this quantile score
			sort(quantiles[i].begin(), quantiles[i].end(), compareQuanMembers);
			
			float high = quantiles[i][int(quantiles[i].size() * 0.99)].score;
			float low =  quantiles[i][int(quantiles[i].size() * 0.01)].score;
			
			vector<quanMember> temp;
			
			//look at each value in quantiles to see if it is an outlier
			for (int j = 0; j < quantiles[i].size(); j++) {
				//is this score between 1 and 99%
				if ((quantiles[i][j].score > low) && (quantiles[i][j].score < high)) {
					temp.push_back(quantiles[i][j]);
				}
			}
			
			quantiles[i] = temp;
		}

/*
		//find contributer with most offending score related to it
		int largestContrib = findLargestContrib(seen);
	
		//while you still have guys to eliminate
		while (contributions.size() > 0) {
		
			mothurOut("Removing scores contributed by sequence " + toString(largestContrib) + " in your template file."); mothurOutEndLine();
			
			//remove from quantiles all scores that were made using this largestContrib
			for (int i = 0; i < quantiles.size(); i++) {
//cout << "before remove " << quantiles[i].size() << '\t';
				removeContrib(largestContrib, quantiles[i]);
//cout << "after remove " << quantiles[i].size() << endl;
			}
//cout << count << " largest contrib = " << largestContrib << endl;  count++;
			//remove from contributions all scores that were made using this largestContrib
			removeContrib(largestContrib, contributions);
			
			//"erase" largestContrib
			seen[largestContrib] = -1;
			
			//get next largestContrib
			largestContrib = findLargestContrib(seen);
		}
ABOVE IS ATTEMPT #1		
**************************************************************************************************
BELOW IS ATTEMPT #2		
		vector<int> marked = returnObviousOutliers(quantiles, num);
		vector<int> copyMarked = marked;
		
		//find the 99% of marked
		sort(copyMarked.begin(), copyMarked.end());
		int high = copyMarked[int(copyMarked.size() * 0.99)];
cout << "high = " << high << endl;		
		
		for(int i = 0; i < marked.size(); i++) {
			if (marked[i] > high) { 
				mothurOut("Removing scores contributed by sequence " + toString(marked[i]) + " in your template file."); mothurOutEndLine();
				for (int i = 0; i < quantiles.size(); i++) {
					removeContrib(marked[i], quantiles[i]);
				}
			}

		}
		
		
		//adjust quantiles
		for (int i = 0; i < quantiles.size(); i++) {
			vector<float> temp;
			
			if (quantiles[i].size() == 0) {
				//in case this is not a distance found in your template files
				for (int g = 0; g < 6; g++) {
					temp.push_back(0.0);
				}
			}else{
				
				sort(quantiles[i].begin(), quantiles[i].end(), compareQuanMembers);
				
				//save 10%
				temp.push_back(quantiles[i][int(quantiles[i].size() * 0.10)].score);
				//save 25%
				temp.push_back(quantiles[i][int(quantiles[i].size() * 0.25)].score);
				//save 50%
				temp.push_back(quantiles[i][int(quantiles[i].size() * 0.5)].score);
				//save 75%
				temp.push_back(quantiles[i][int(quantiles[i].size() * 0.75)].score);
				//save 95%
				temp.push_back(quantiles[i][int(quantiles[i].size() * 0.95)].score);
				//save 99%
				temp.push_back(quantiles[i][int(quantiles[i].size() * 0.99)].score);
				
			}
			
			quan[i] = temp;
			
		}
*/
		
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "removeObviousOutliers");
		exit(1);
	}
}
//***************************************************************************************************************
//put quanMember in the vector based on how far they are from the 99% or 1%.  Biggest offenders in front.
/*vector<quanMember> DeCalculator::sortContrib(map<quanMember*, float> quan) {
	try{
		
		vector<quanMember> newQuan;
		map<quanMember*, float>::iterator it;
		
		while (quan.size() > 0) {
			
			 map<quanMember*, float>::iterator largest = quan.begin(); 
			  
			//find biggest member
			for (it = quan.begin(); it != quan.end(); it++) {
				if (it->second > largest->second) {  largest = it;  }
			}
cout << largest->second << '\t' << largest->first->score << '\t' << largest->first->member1 << '\t' << largest->first->member2 << endl;
			//save it 
			newQuan.push_back(*(largest->first));
		
			//erase from quan
			quan.erase(largest);
		}
		
		return newQuan;
		
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "sortContrib");
		exit(1);
	}
}

//***************************************************************************************************************
//used by removeObviousOutliers which was attempt to increase sensitivity of chimera detection...not currently used...
int DeCalculator::findLargestContrib(vector<int> seen) {
	try{
		
		int largest = 0;
		
		int largestContribs;
		
		for (int i = 0; i < seen.size(); i++)  {  
			
			if (seen[i] > largest) {
				largestContribs = i;
				largest = seen[i];
			}
		}
		
		return largestContribs;
		
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "findLargestContribs");
		exit(1);
	}
}
//***************************************************************************************************************
void DeCalculator::removeContrib(int bad, vector<quanMember>& quan) {
	try{
	
		vector<quanMember> newQuan;
		for (int i = 0; i < quan.size(); i++)  {  
			if ((quan[i].member1 != bad) && (quan[i].member2 != bad) ) {  
				newQuan.push_back(quan[i]);  
			}
		}
		
		quan = newQuan;
		
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "removeContrib");
		exit(1);
	}
}
*/
//***************************************************************************************************************
float DeCalculator::findAverage(vector<float> myVector) {
	try{
		
		float total = 0.0;
		for (int i = 0; i < myVector.size(); i++)  {  total += myVector[i];  }
		
		float average = total / (float) myVector.size();
		
		return average;
		
	}
	catch(exception& e) {
		errorOut(e, "DeCalculator", "findAverage");
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
		errorOut(e, "DeCalculator", "getCoef");
		exit(1);
	}
}
//***************************************************************************************************************


