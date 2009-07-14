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

Pintail::Pintail(string name) {
	try {
		fastafile = name;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "Pintail");
		exit(1);
	}
}
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
		
		for (itCoef = DE.begin(); itCoef != DE.end(); itCoef++) {
			
			out << itCoef->first->getName() << '\t' << itCoef->second << endl;
			out << "Observed\t";
			
			itObsDist = obsDistance.find(itCoef->first);
			for (int i = 0; i < itObsDist->second.size(); i++) {  out << itObsDist->second[i] << '\t';  }
			out << endl;
			
			out << "Expected\t";
			
			itExpDist = expectedDistance.find(itCoef->first);
			for (int i = 0; i < itExpDist->second.size(); i++) {  out << itExpDist->second[i] << '\t';  }
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
		
		distCalculator = new ignoreGaps();
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		templateSeqs = readSeqs(templateFile);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		//if window is set to default
		if (window == 0) {  if (querySeqs[0]->getAligned().length() > 800) {  setWindow(200); }
							else{  setWindow((querySeqs[0]->getAligned().length() / 4)); }  } 
		else if (window > (querySeqs[0]->getAligned().length() / 4)) { 
				mothurOut("You have selected to large a window size for you sequences.  I will choose a smaller window."); mothurOutEndLine();
				setWindow((querySeqs[0]->getAligned().length() / 4)); 
		}
	
		//calculate number of iters 
		iters = (querySeqs[0]->getAligned().length() - window + 1) / increment;
cout << "length = " << querySeqs[0]->getAligned().length() << " window = " << window << " increment = " << increment << " iters = " << iters << endl;			
		int linesPerProcess = processors / numSeqs;
		
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
		
		//map query sequences to their most similiar sequences in the template - Parallelized
		mothurOut("Finding closest sequence in template to each sequence... "); cout.flush();
		if (processors == 1) {   findPairs(lines[0]->start, lines[0]->end);  } 
		else {		createProcessesPairs();		}
		mothurOut("Done."); mothurOutEndLine();
		
		//find Oqs for each sequence - the observed distance in each window - Parallelized
		mothurOut("Calculating observed percentage differences for each sequence... "); cout.flush();
		if (processors == 1) {   calcObserved(lines[0]->start, lines[0]->end);  } 
		else {		createProcessesObserved();		}
		mothurOut("Done."); mothurOutEndLine();
		
		//find P
		mothurOut("Calculating expected percentage differences for each sequence... "); cout.flush();
		vector<float> probabilityProfile = calcFreq(templateSeqs);
		
		//make P into Q
		for (int i = 0; i < probabilityProfile.size(); i++)  {	probabilityProfile[i] = 1 - probabilityProfile[i];	}

		//find Qav
		averageProbability = findQav(probabilityProfile);
	
		//find Coefficient - maps a sequence to its coefficient
		seqCoef = getCoef(averageProbability);
	
		//find Eqs for each sequence - the expected distance in each window - Parallelized
		if (processors == 1) {   calcExpected(lines[0]->start, lines[0]->end);  } 
		else {		createProcessesExpected();		}
		mothurOut("Done."); mothurOutEndLine();
		
		//find deviation - Parallelized
		mothurOut("Finding deviation from expected... "); cout.flush();
		if (processors == 1) {   calcDE(lines[0]->start, lines[0]->end);  } 
		else {		createProcessesDE();		}
		mothurOut("Done."); mothurOutEndLine();
		
				
		//free memory
		for (int i = 0; i < lines.size(); i++)			{	delete lines[i];		}
		delete distCalculator;	

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
			
			if (current->getAligned() == "") { current->setAligned(current->getUnaligned()); }
			//takes out stuff is needed
			current->setUnaligned(current->getUnaligned());
			
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
//calculate the distances from each query sequence to all sequences in the template to find the closest sequence
void Pintail::findPairs(int start, int end) {
	try {
		
		for(int i = start; i < end; i++){
		
			float smallest = 10000.0;
			Sequence query = *(querySeqs[i]);
		
			for(int j = 0; j < templateSeqs.size(); j++){
				
				Sequence temp = *(templateSeqs[j]);
				
				distCalculator->calcDist(query, temp);
				float dist = distCalculator->getDist();
				
				if (dist < smallest) { 
				
					bestfit[querySeqs[i]] = templateSeqs[j];
					smallest = dist;
				}
			}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "findPairs");
		exit(1);
	}
}
//***************************************************************************************************************
void Pintail::calcObserved(int start, int end) {
	try {
	
						
		for(int i = start; i < end; i++){
		
			itBest = bestfit.find(querySeqs[i]);
			Sequence* query;
			Sequence* subject;
		
			if (itBest != bestfit.end()) {
				query = itBest->first;
				subject = itBest->second;
			}else{ mothurOut("Error in calcObserved"); mothurOutEndLine(); }
//cout << query->getName() << '\t' << subject->getName() << endl;			
			
			int startpoint = 0; 
			for (int m = 0; m < iters; m++) {

				string seqFrag = query->getAligned().substr(startpoint, window);
				string seqFragsub = subject->getAligned().substr(startpoint, window);
								
				int diff = 0;
                for (int b = 0; b < seqFrag.length(); b++) {
                  
                    //if this is not a gap
                    if ((isalpha(seqFrag[b])) && (isalpha(seqFragsub[b]))) {
                        //and they are different - penalize
                        if (seqFrag[b] != seqFragsub[b]) { diff++; }
                    }
                }
               
                //percentage of mismatched bases
                float dist = diff / (float)seqFrag.length();       
				
				obsDistance[query].push_back(dist);
				
				startpoint += increment;
			}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcObserved");
		exit(1);
	}
}

//***************************************************************************************************************
void Pintail::calcExpected(int start, int end) {
	try {
		
	
		//for each sequence
		for(int i = start; i < end; i++){
			
			itCoef = seqCoef.find(querySeqs[i]);
			float coef = itCoef->second;
			
			//for each window
			vector<float> queryExpected;
			for (int m = 0; m < iters; m++) {		
				float expected = averageProbability[m] * coef;
				queryExpected.push_back(expected);	
//cout << "average variabilty over window = " << averageProbability[m] << " coef = " << coef << " ei = "  << expected << '\t' <<  "window = " << m << endl;
			}
			
			expectedDistance[querySeqs[i]] = queryExpected;
			
		}
				
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcExpected");
		exit(1);
	}
}
//***************************************************************************************************************
void Pintail::calcDE(int start, int end) {
	try {
		
	
		//for each sequence
		for(int i = start; i < end; i++){
			
			itObsDist = obsDistance.find(querySeqs[i]);
			vector<float> obs = itObsDist->second;
			
			itExpDist = expectedDistance.find(querySeqs[i]);
			vector<float> exp = itExpDist->second;
//	cout << "difference between obs and exp = " << abs(obs[m] - exp[m]) << endl;	
			//for each window
			float sum = 0.0;  //sum = sum from 1 to m of (oi-ei)^2
			for (int m = 0; m < iters; m++) { 		sum += ((obs[m] - exp[m]) * (obs[m] - exp[m]));		}
			
			float de = sqrt((sum / (iters - 1)));
			
			DE[querySeqs[i]] = de;
		}
				
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
			for (int m = 0; m < freq.size(); m++) {    if (freq[m] > highest) {  highest = freq[m];  }		}
			
			//add in gaps - so you can effectively "ignore them"
			highest += gaps;
			
			float highFreq = highest / (float) seqs.size();	
			
			float Pi;
			Pi =  (highFreq - 0.25) / 0.75;  
				
			prob.push_back(Pi); 
		}
		
		return prob;
				
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcFreq");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float> Pintail::findQav(vector<float> prob) {
	try {
		vector<float> averages;
		
		//for each window find average
		int startpoint = 0;
		for (int m = 0; m < iters; m++) {
			
			float average = 0.0;
			for (int i = startpoint; i < (startpoint+window); i++) {   average += prob[i];	}
			
			average = average / window;
//cout << average << endl;			
			//save this windows average
			averages.push_back(average);
		
			startpoint += increment;
		}
		
		return averages;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "findQav");
		exit(1);
	}
}
//***************************************************************************************************************
map<Sequence*, float> Pintail::getCoef(vector<float> prob) {
	try {
		map<Sequence*, float> coefs;
		
		//find average prob
		float probAverage = 0.0;
		for (int i = 0; i < prob.size(); i++) {   probAverage += prob[i];	}
		probAverage = probAverage / (float) prob.size();
cout << "(sum of ai) / m = " << probAverage << endl;		
		//find a coef for each sequence
		map<Sequence*, vector<float> >::iterator it;
		for (it = obsDistance.begin(); it != obsDistance.end(); it++) {
			
			vector<float> temp = it->second;
			Sequence* tempSeq = it->first;
			
			//find observed average 
			float obsAverage = 0.0;
			for (int i = 0; i < temp.size(); i++) {   obsAverage += temp[i];	}
			obsAverage = obsAverage / (float) temp.size();
cout << tempSeq->getName() << '\t' << obsAverage << endl;			
			float coef = obsAverage / probAverage;
cout  << tempSeq->getName() << '\t' << "coef = " << coef << endl;			
			//save this sequences coefficient
			coefs[tempSeq] = coef;
		}
		
						
		return coefs;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "getCoef");
		exit(1);
	}
}

/**************************************************************************************************/

void Pintail::createProcessesPairs() {
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
				findPairs(lines[process]->start, lines[process]->end);
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
#else
		findPairs(lines[0]->start, lines[0]->end);

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcessesPairs");
		exit(1);
	}
}

/**************************************************************************************************/

void Pintail::createProcessesObserved() {
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
				calcObserved(lines[process]->start, lines[process]->end);
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
#else
		calcObserved(lines[0]->start, lines[0]->end);

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcessesObserved");
		exit(1);
	}
}

//***************************************************************************************************************

void Pintail::createProcessesExpected() {
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
				calcExpected(lines[process]->start, lines[process]->end);
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
#else
		calcExpected(lines[0]->start, lines[0]->end);

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcessesExpected");
		exit(1);
	}
}

/**************************************************************************************************/

void Pintail::createProcessesDE() {
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
				calcDE(lines[process]->start, lines[process]->end);
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
#else
		calcDE(lines[0]->start, lines[0]->end);

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcessesDE");
		exit(1);
	}
}

//***************************************************************************************************************


