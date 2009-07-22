/*
 *  pintail.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "pintail.h"
#include "ignoregaps.h"

//***************************************************************************************************************

Pintail::Pintail(string filename, string temp) {  fastafile = filename;  templateFile = temp;  }
//***************************************************************************************************************

Pintail::~Pintail() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)		{  delete querySeqs[i];		}
		for (int i = 0; i < templateSeqs.size(); i++)	{  delete templateSeqs[i];	}
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "~Pintail");
		exit(1);
	}
}	
//***************************************************************************************************************
void Pintail::print(ostream& out) {
	try {
		
		for (int i = 0; i < querySeqs.size(); i++) {
			
			out << querySeqs[i]->getName() << '\t' << "div: " << deviation[i] << "\tstDev: " << DE[i] << endl;
			out << "Observed\t";
			
			for (int j = 0; j < obsDistance[i].size(); j++) {  out << obsDistance[i][j] << '\t';  }
			out << endl;
			
			out << "Expected\t";
			
			for (int m = 0; m < expectedDistance[i].size(); m++) {  out << expectedDistance[i][m] << '\t';  }
			out << endl;
			
		}
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "print");
		exit(1);
	}
}

//***************************************************************************************************************
void Pintail::getChimeras() {
	try {
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		templateSeqs = readSeqs(templateFile);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		obsDistance.resize(numSeqs);
		expectedDistance.resize(numSeqs);
		seqCoef.resize(numSeqs);
		DE.resize(numSeqs);
		Qav.resize(numSeqs);
		bestfit.resize(numSeqs);
		deviation.resize(numSeqs);
		trimmed.resize(numSeqs);
		windowSizes.resize(numSeqs, window);
		windows.resize(numSeqs);
		
		//break up file if needed
		int linesPerProcess = processors / numSeqs;
		
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
		
		distcalculator = new ignoreGaps();

				
		if (processors == 1) { 
			mothurOut("Finding closest sequence in template to each sequence... "); cout.flush();
			bestfit = findPairs(lines[0]->start, lines[0]->end);
/*for (int m = 0; m < templateSeqs.size(); m++)  {
	if (templateSeqs[m]->getName() == "159481") {  bestfit[17] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "100137") {  bestfit[16] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "112956") {  bestfit[15] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "102326") {  bestfit[14] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "66229") {  bestfit[13] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "206276") {  bestfit[12] = *(templateSeqs[m]); }
    if (templateSeqs[m]->getName() == "63607") {  bestfit[11] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "7056") {  bestfit[10] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "7088") {  bestfit[9] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "17553") {  bestfit[8] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "131723") {  bestfit[7] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "69013") {  bestfit[6] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "24543") {  bestfit[5] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "27824") {  bestfit[4] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "1456") {  bestfit[3] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "1456") {  bestfit[2] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "141312") {  bestfit[1] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "141312") {  bestfit[0] = *(templateSeqs[m]); }


}*/
			
			for (int j = 0; j < bestfit.size(); j++) { 
				//chops off beginning and end of sequences so they both start and end with a base
				trimSeqs(querySeqs[j], bestfit[j], j);  
			}
			mothurOut("Done."); mothurOutEndLine();
			
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				it = trimmed[i].begin();
				vector<int> win = findWindows(querySeqs[i], it->first, it->second, windowSizes[i]);
				windows[i] = win;
			}

			
		} else {		createProcessesSpots();		}

		//find P
		if (consfile == "") {   probabilityProfile = calcFreq(templateSeqs);  }
		else				{   probabilityProfile = readFreq();			  }
		
		//make P into Q
		for (int i = 0; i < probabilityProfile.size(); i++)  {	probabilityProfile[i] = 1 - probabilityProfile[i];	}
		
		if (processors == 1) { 
						
			mothurOut("Calculating observed distance... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				vector<float> obsi = calcObserved(querySeqs[i], bestfit[i], windows[i], windowSizes[i]);
				obsDistance[i] = obsi;
			}
			mothurOut("Done."); mothurOutEndLine();
			
			
			
			mothurOut("Finding variability... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				vector<float> q = findQav(windows[i], windowSizes[i]);
				Qav[i] = q;
			}
			mothurOut("Done."); mothurOutEndLine();
			
			
			
			mothurOut("Calculating alpha... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				float alpha = getCoef(obsDistance[i], Qav[i]);
				seqCoef.push_back(alpha);
			}
			mothurOut("Done."); mothurOutEndLine();
		
		
		
			mothurOut("Calculating expected distance... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				vector<float> exp = calcExpected(Qav[i], seqCoef[i]);
				expectedDistance[i] = exp;
			}
			mothurOut("Done."); mothurOutEndLine();
			
			
			
			mothurOut("Finding deviation... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				float de = calcDE(obsDistance[i], expectedDistance[i]);
				DE[i] = de;
				
				it = trimmed[i].begin();
				float dist = calcDist(querySeqs[i], bestfit[i], it->first, it->second); 
				deviation[i] = dist;
			}
			mothurOut("Done."); mothurOutEndLine();
			
		} 
		else {		createProcesses();		}
		
		delete distcalculator;
		
		//quantiles are used to determine whether the de values found indicate a chimera
		//if you have to calculate them, its time intensive because you are finding the de and deviation values for each 
		//combination of sequences in the template
		if (quanfile != "") {   readQuantiles();  }
		else {
			if (processors == 1) { 
				quantiles = getQuantiles(lines[0]->start, lines[0]->end);
			}else {		createProcessesQuan();		}
		}

	
		//free memory
		for (int i = 0; i < lines.size(); i++)			{	delete lines[i];		}
			

	}
	catch(exception& e) {
		errorOut(e, "Pintail", "getChimeras");
		exit(1);
	}
}

//***************************************************************************************************************

vector<Sequence*> Pintail::readSeqs(string file) {
	try {
	
		ifstream in;
		openInputFile(file, in);
		vector<Sequence*> container;
		
		//read in seqs and store in vector
		while(!in.eof()){
			Sequence* current = new Sequence(in);
			
			//if (current->getAligned() == "") { current->setAligned(current->getUnaligned()); }
			//takes out stuff is needed
			//current->setUnaligned(current->getUnaligned());
			
			container.push_back(current);
			
			gobble(in);
		}
		
		in.close();
		return container;
		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "readSeqs");
		exit(1);
	}
}

//***************************************************************************************************************
//num is query's spot in querySeqs
map<int, int> Pintail::trimSeqs(Sequence* query, Sequence subject, int num) {
	try {
		
		string q = query->getAligned();
		string s = subject.getAligned();
		
		int front = 0;
		for (int i = 0; i < q.length(); i++) {
			if (isalpha(q[i]) && isalpha(s[i])) { front = i; break;  }
		}
		
		int back = 0;		
		for (int i = q.length(); i >= 0; i--) {
			if (isalpha(q[i]) && isalpha(s[i])) { back = i; break;  }
		}
		
		//if num = -1 then you are calling this from quantiles
		if (num != -1) { 
			trimmed[num][front] = back; 
			return trimmed[num];
		}
		
		map<int, int> temp;
		temp[front] = back;
		return temp;
		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "trimSeqs");
		exit(1);
	}
}

//***************************************************************************************************************

vector<float> Pintail::readFreq() {
	try {
	
		ifstream in;
		openInputFile(consfile, in);
		
		vector<float> prob;
		
		//read in probabilities and store in vector
		int pos; float num; 
		
		while(!in.eof()){
			
			in >> pos >> num;
			
			//do you want this spot
			prob.push_back(num);  
			
			gobble(in);
		}
		
		in.close();
		return prob;
		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "readFreq");
		exit(1);
	}
}

//***************************************************************************************************************

vector< vector<float> > Pintail::readQuantiles() {
	try {
	
		ifstream in;
		openInputFile(quanfile, in);
		
		vector< vector<float> > quan;
		
		//read in probabilities and store in vector
		int num; float ten, twentyfive, fifty, seventyfive, ninetyfive, ninetynine; 
		
		while(!in.eof()){
			
			in >> num >> ten >> twentyfive >> fifty >> seventyfive >> ninetyfive >> ninetynine; 
			
			vector <float> temp;
			
			temp.push_back(ten); 
			temp.push_back(twentyfive);
			temp.push_back(fifty);
			temp.push_back(seventyfive);
			temp.push_back(ninetyfive);
			temp.push_back(ninetynine);
			
			//do you want this spot
			quan.push_back(temp);  
			
			gobble(in);
		}
		
		in.close();
		return quan;
		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "readQuantiles");
		exit(1);
	}
}


//***************************************************************************************************************
//calculate the distances from each query sequence to all sequences in the template to find the closest sequence
vector<Sequence> Pintail::findPairs(int start, int end) {
	try {
		
		vector<Sequence> seqsMatches;  seqsMatches.resize(end-start);
		
		for(int i = start; i < end; i++){
		
			float smallest = 10000.0;
			Sequence query = *(querySeqs[i]);
		
			for(int j = 0; j < templateSeqs.size(); j++){
				
				Sequence temp = *(templateSeqs[j]);
				
				distcalculator->calcDist(query, temp);
				float dist = distcalculator->getDist();
				
				if (dist < smallest) { 
					seqsMatches[i] = *(templateSeqs[j]);
					smallest = dist;
				}
			}
		}
		
		return seqsMatches;
	
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "findPairs");
		exit(1);
	}
}

//***************************************************************************************************************
//find the window breaks for each sequence - this is so you can move ahead by bases.
vector<int>  Pintail::findWindows(Sequence* query, int front, int back, int& size) {
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
		errorOut(e, "Pintail", "findWindows");
		exit(1);
	}
}

//***************************************************************************************************************
vector<float> Pintail::calcObserved(Sequence* query, Sequence subject, vector<int> window, int size) {
	try {
		
		vector<float> temp;
				
		int startpoint = 0; 
		for (int m = 0; m < windows.size(); m++) {
						
			string seqFrag = query->getAligned().substr(window[startpoint], size);
			string seqFragsub = subject.getAligned().substr(window[startpoint], size);
								
			int diff = 0;
			for (int b = 0; b < seqFrag.length(); b++) {
				if (seqFrag[b] != seqFragsub[b]) { diff++; }
			}
               
			//percentage of mismatched bases
			float dist;
			dist = diff / (float) seqFrag.length() * 100;       
				
			temp.push_back(dist);
				
			startpoint++;
		}
			
		return temp;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcObserved");
		exit(1);
	}
}
//***************************************************************************************************************
float Pintail::calcDist(Sequence* query, Sequence subject, int front, int back) {
	try {
		
		//so you only look at the trimmed part of the sequence
		int cutoff = back - front;  
			
		//from first startpoint with length back-front
		string seqFrag = query->getAligned().substr(front, cutoff);
		string seqFragsub = subject.getAligned().substr(front, cutoff);
														
		int diff = 0;
		for (int b = 0; b < seqFrag.length(); b++) {
			if (seqFrag[b] != seqFragsub[b]) { diff++; }
		}
               
		//percentage of mismatched bases
		float dist = diff / (float) seqFrag.length() * 100;       
				
		return dist;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcDist");
		exit(1);
	}
}

//***************************************************************************************************************
vector<float> Pintail::calcExpected(vector<float> qav, float coef) {
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
		errorOut(e, "Pintail", "calcExpected");
		exit(1);
	}
}
//***************************************************************************************************************
float Pintail::calcDE(vector<float> obs, vector<float> exp) {
	try {
		
		//for each window
		float sum = 0.0;  //sum = sum from 1 to m of (oi-ei)^2
		for (int m = 0; m < obsDistance.size(); m++) { 		sum += ((obs[m] - exp[m]) * (obs[m] - exp[m]));		}
			
		float de = sqrt((sum / (obsDistance.size() - 1)));
			
		return de;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcDE");
		exit(1);
	}
}

//***************************************************************************************************************

vector<float> Pintail::calcFreq(vector<Sequence*> seqs) {
	try {

		vector<float> prob;
		string freqfile = getRootName(templateFile) + "probability";
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
			if ( (seqs.size() - gaps) == 0 ) {  highFreq = 1.0;  }			
			else { highFreq = highest / (float) (seqs.size() - gaps);	 }
			//highFreq = highest / (float) seqs.size();
			//cout << i << '\t' << highFreq << endl;
			
			float Pi;
			Pi =  (highFreq - 0.25) / 0.75; 
			
			//cannot have probability less than 0.
			if (Pi < 0) { Pi = 0.0; }
			
			//saves this for later
			outFreq << i+1 << '\t' << Pi << endl;
			
			prob.push_back(Pi); 
		}
		
		outFreq.close();
		
		return prob;
				
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcFreq");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float>  Pintail::findQav(vector<int> window, int size) {
	try {
		vector<float>  averages; 
				
		//for each window find average
		for (int m = 0; m < windows.size(); m++) {
				
			float average = 0.0;
				
			//while you are in the window for this sequence
			for (int j = window[m]; j < (window[m]+size); j++) {   average += probabilityProfile[j];	}
				
			average = average / size;
	
			//save this windows average
			averages.push_back(average);
		}
				
		return averages;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "findQav");
		exit(1);
	}
}

//***************************************************************************************************************
vector< vector<float> > Pintail::getQuantiles(int start, int end) {
	try {
		vector< vector<float> > quan; 
		
		//for each sequence
		for(int i = start; i < end; i++){
		
			Sequence* query = templateSeqs[i];
			
			//compare to every other sequence in template
			for(int j = 0; j < i; j++){
				
				Sequence subject = *(templateSeqs[j]);
				
				map<int, int> trim = trimSeqs(query, subject, -1);
				
				
				
			
			
			}
			
			
			
		}
		return quan;
						
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "findQav");
		exit(1);
	}
}

//***************************************************************************************************************
float Pintail::getCoef(vector<float> obs, vector<float> qav) {
	try {
	
		//find average prob for this seqs windows
		float probAverage = 0.0;
		for (int j = 0; j < qav.size(); j++) {   probAverage += qav[j];	}
		probAverage = probAverage / (float) Qav.size();
		
		//find observed average 
		float obsAverage = 0.0;
		for (int j = 0; j < obs.size(); j++) {   obsAverage += obs[j];	}
		obsAverage = obsAverage / (float) obs.size();
		

		float coef = obsAverage / probAverage;
						
		return coef;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "getCoef");
		exit(1);
	}
}
/**************************************************************************************************/

void Pintail::createProcessesSpots() {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		vector<int> processIDS;
		vector< vector<int> > win; win.resize(querySeqs.size());
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
				
				vector<Sequence> tempbest;
				tempbest = findPairs(lines[process]->start, lines[process]->end);
				int count = 0;
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					bestfit[i] = tempbest[count];
					
					//chops off beginning and end of sequences so they both start and end with a base
					trimSeqs(querySeqs[i], bestfit[i], i);
					count++;
				}
				
				
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					vector<int> temp = findWindows(querySeqs[i], it->first, it->second, windowSizes[i]);
					win[i] = temp;
				}
				
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		windows = win;
#else
		bestfit = findPairs(lines[0]->start, lines[0]->end);
		for (int j = 0; j < bestfit.size(); j++) {
				//chops off beginning and end of sequences so they both start and end with a base
				trimSeqs(querySeqs[j], bestfit[j], j);  
		}

		for (int i = lines[0]->start; i < lines[0]->end; i++) {
				it = trimmed[i].begin();
				map<int, int> win = findWindows(querySeqs[i], it->first, it->second, windowSizes[i]);
				windows[i] = win;
		}

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcessesSpots");
		exit(1);
	}
}


/**************************************************************************************************/

void Pintail::createProcesses() {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		vector<int> processIDS;
		
		vector< vector<float> > exp;  exp.resize(querySeqs.size());
		vector<float> de; de.resize(querySeqs.size());
		vector< vector<float> > obs; obs.resize(querySeqs.size());
		vector<float> dev; dev.resize(querySeqs.size());
		
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
				
				//calc obs
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					vector<float> obsi = calcObserved(querySeqs[i], bestfit[i], windows[i], windowSizes[i]);
					obs[i] = obsi;
				
					//calc Qav
					vector<float> q = findQav(windows[i], windowSizes[i]);
					
					//get alpha
					float alpha = getCoef(obsDistance[i], q);
					
					//find expected
					vector<float> exp = calcExpected(q, alpha);
					expectedDistance[i] = exp;
					
					//get de and deviation
					float dei = calcDE(obsi, exp);
					de[i] = dei;
					
					it = trimmed[i].begin();
					float dist = calcDist(querySeqs[i], bestfit[i], it->first, it->second); 
					dev[i] = dist;
				}
				
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		obsDistance = obs;
		expectedDistance = exp;
		DE = de;
		deviation = dev;
		
#else
			mothurOut("Calculating observed distance... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				vector<float> obsi = calcObserved(querySeqs[i], bestfit[i], windows[i], windowSizes[i]);
				obsDistance[i] = obsi;
			}
			mothurOut("Done."); mothurOutEndLine();
			
			
			
			mothurOut("Finding variability... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				vector<float> q = findQav(windows[i], windowSizes[i]);
				Qav[i] = q;
			}
			mothurOut("Done."); mothurOutEndLine();
			
			
			
			mothurOut("Calculating alpha... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				float alpha = getCoef(obsDistance[i], Qav[i]);
				seqCoef.push_back(alpha);
			}
			mothurOut("Done."); mothurOutEndLine();
		
		
		
			mothurOut("Calculating expected distance... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				vector<float> exp = calcExpected(Qav[i], seqCoef[i]);
				expectedDistance[i] = exp;
			}
			mothurOut("Done."); mothurOutEndLine();
			
			
			
			mothurOut("Finding deviation... "); cout.flush();
			for (int i = lines[0]->start; i < lines[0]->end; i++) {
				float de = calcDE(obsDistance[i], expectedDistance[i]);
				DE[i] = de;
				
				it = trimmed[i].begin();
				float dist = calcDist(querySeqs[i], bestfit[i], it->first, it->second); 
				deviation[i] = dist;
			}
			mothurOut("Done."); mothurOutEndLine();

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcesses");
		exit(1);
	}
}


/**************************************************************************************************/

void Pintail::createProcessesQuan() {
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
				
				
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
#else

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcessesQuan");
		exit(1);
	}
}


//***************************************************************************************************************


