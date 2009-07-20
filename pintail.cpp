/*
 *  pintail.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "pintail.h"
#include "eachgapdist.h"

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
		trim.resize(numSeqs);
		deviation.resize(numSeqs);
		windowSizes.resize(numSeqs, window);
		
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
		
		distcalculator = new eachGapDist();
		
		if (processors == 1) { 
			mothurOut("Finding closest sequence in template to each sequence... "); cout.flush();
			bestfit = findPairs(lines[0]->start, lines[0]->end);
for (int m = 0; m < templateSeqs.size(); m++)  {
	if (templateSeqs[m]->getName() == "198806") {  bestfit[0] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "198806") {  bestfit[1] = *(templateSeqs[m]); }
	if (templateSeqs[m]->getName() == "108139") {  bestfit[2] = *(templateSeqs[m]); }
}
			
for (int j = 0; j < bestfit.size(); j++) {//cout << querySeqs[j]->getName() << '\t' << "length = " <<  querySeqs[j]->getAligned().length() << '\t' << bestfit[j].getName() << " length = " <<  bestfit[j].getAligned().length() <<  endl; 
				//chops off beginning and end of sequences so they both start and end with a base
				trimSeqs(querySeqs[j], bestfit[j], j);  
//cout << "NEW SEQ PAIR" << querySeqs[j]->getAligned() << endl << "IN THE MIDDLE" <<  endl << bestfit[j].getAligned() << endl; 

}

			mothurOut("Done."); mothurOutEndLine();

			windows = findWindows(lines[0]->start, lines[0]->end);
		} else {		createProcessesSpots();		}

		//find P
		if (consfile == "") {   probabilityProfile = calcFreq(templateSeqs);  }
		else				{   probabilityProfile = readFreq();			  }
		
		//make P into Q
		for (int i = 0; i < probabilityProfile.size(); i++)  {	probabilityProfile[i] = 1 - probabilityProfile[i];	}
		
		if (processors == 1) { 
						
			mothurOut("Calculating observed distance... "); cout.flush();
			obsDistance = calcObserved(lines[0]->start, lines[0]->end);
			mothurOut("Done."); mothurOutEndLine();
			
			mothurOut("Finding variability... "); cout.flush();
			Qav = findQav(lines[0]->start, lines[0]->end);
for (int i = 0; i < Qav.size(); i++) {
cout << querySeqs[i]->getName() << " = ";
for (int u = 0; u < Qav[i].size();u++) {   cout << Qav[i][u] << '\t';  }
cout << endl << endl;
}


			mothurOut("Done."); mothurOutEndLine();
			
			mothurOut("Calculating alpha... "); cout.flush();
			seqCoef = getCoef(lines[0]->start, lines[0]->end);
for (int i = 0; i < seqCoef.size(); i++) {
cout << querySeqs[i]->getName() << " coef = " << seqCoef[i] << endl;
}

			mothurOut("Done."); mothurOutEndLine();
			
			mothurOut("Calculating expected distance... "); cout.flush();
			expectedDistance = calcExpected(lines[0]->start, lines[0]->end);
			mothurOut("Done."); mothurOutEndLine();
			
			mothurOut("Finding deviation... "); cout.flush();
			DE = calcDE(lines[0]->start, lines[0]->end); 
			deviation = calcDist(lines[0]->start, lines[0]->end); 
			mothurOut("Done."); mothurOutEndLine();
			
			
			
		} 
		else {		createProcesses();		}
		
		delete distcalculator;
	
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
//num is query's spot in querySeqs
void Pintail::trimSeqs(Sequence* query, Sequence& subject, int num) {
	try {
	
		string q = query->getAligned();
		string s = subject.getAligned();
	
		int front = 0;
		for (int i = 0; i < q.length(); i++) {
			if (isalpha(q[i]) && isalpha(s[i])) { front = i; break;  }
		}
		
		q = q.substr(front, q.length());
		s = s.substr(front, s.length());
		
		int back = 0;		
		for (int i = q.length(); i >= 0; i--) {
			if (isalpha(q[i]) && isalpha(s[i])) { back = i; break;  }
		}
	
		q = q.substr(0, back);
		s = s.substr(0, back);

		trim[num][front] = back;
	
		//save 
		query->setAligned(q);
		query->setUnaligned(q);
		subject.setAligned(s);
		subject.setUnaligned(s);
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
//find the window breaks for each sequence
vector< vector<int> > Pintail::findWindows(int start, int end) {
	try {
		
		vector< vector<int> > win;  win.resize(end-start);
		
		//for each sequence
		int count = 0;
		for(int i = start; i < end; i++){
			
			//if window is set to default
			if (windowSizes[i] == 0) {  if (querySeqs[i]->getAligned().length() > 1200) {  windowSizes[i] = 300; }
							else{  windowSizes[i] = (querySeqs[i]->getAligned().length() / 4); }  } 
			else if (windowSizes[i] > (querySeqs[i]->getAligned().length() / 4)) { 
					mothurOut("You have selected to large a window size for sequence " + querySeqs[i]->getName() + ".  I will choose an appropriate window size."); mothurOutEndLine();
					windowSizes[i] = (querySeqs[i]->getAligned().length() / 4); 
			}
	
	//cout << "length = " << querySeqs[i]->getAligned().length() << " window = " << windowSizes[i] << " increment = " << increment << endl;			
				

			string seq = querySeqs[i]->getAligned();
			int numBases = querySeqs[i]->getUnaligned().length();
			int spot = 0;
			
			//find location of first base
			for (int j = 0; j < seq.length(); j++) {
				if (isalpha(seq[j])) { spot = j;  break;  }
			}
			
			//save start of seq
			win[count].push_back(spot);
			
			
			//move ahead increment bases at a time until all bases are in a window
			int countBases = 0;
			int totalBases = 0;  //used to eliminate window of blanks at end of sequence
			for (int m = spot; m < seq.length(); m++) {
				
				//count number of bases you see
				if (isalpha(seq[m])) { countBases++; totalBases++;  }
				
				//if you have seen enough bases to make a new window
				if (countBases >= increment) {
					win[count].push_back(m);  //save spot in alignment
					countBases = 0;				//reset bases you've seen in this window
				}
				
				//no need to continue if all your bases are in a window
				if (totalBases == numBases) {   break;  }
			}
			
			count++;
		}
		
		
		
		return win;
	
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "findWindows");
		exit(1);
	}
}

//***************************************************************************************************************
vector< vector<float> > Pintail::calcObserved(int start, int end) {
	try {
		
		vector< vector<float> > temp;
		temp.resize(end-start);
		
		int count = 0;
		for(int i = start; i < end; i++){
		
			Sequence* query = querySeqs[i];
			Sequence subject = bestfit[i];
		
			int startpoint = 0; 
			for (int m = 0; m < windows[i].size(); m++) {

				string seqFrag = query->getAligned().substr(windows[i][startpoint], windowSizes[i]);
				string seqFragsub = subject.getAligned().substr(windows[i][startpoint], windowSizes[i]);
								
				int diff = 0;
                for (int b = 0; b < seqFrag.length(); b++) {
                  
                    //if either the query or subject is not a gap 
                    if ((isalpha(seqFrag[b])) || (isalpha(seqFragsub[b]))) {
                        //and they are different - penalize
                        if (seqFrag[b] != seqFragsub[b]) { diff++; }
                    }
                }
               
                //percentage of mismatched bases
				float dist;
                dist = diff / (float) seqFrag.length() * 100;       
				
				temp[count].push_back(dist);
				
				startpoint++;
			}
			
			count++;
		}
		
		return temp;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcObserved");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float>  Pintail::calcDist(int start, int end) {
	try {
		
		vector<float> temp;
				
		for(int i = start; i < end; i++){
		
			Sequence* query = querySeqs[i];
			Sequence subject = bestfit[i];
			
			string seqFrag = query->getAligned();
			string seqFragsub = subject.getAligned();
														
			int diff = 0;
			for (int b = 0; b < seqFrag.length(); b++) {
                  
				//if either the query or subject is not a gap 
				if ((isalpha(seqFrag[b])) || (isalpha(seqFragsub[b]))) {
					//and they are different - penalize
					if (seqFrag[b] != seqFragsub[b]) { diff++; }
				}
			}
               
			//percentage of mismatched bases
			float dist;
			dist = diff / (float) seqFrag.length() * 100;       
				
			temp.push_back(dist);
		}
		
		return temp;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcDist");
		exit(1);
	}
}

//***************************************************************************************************************
vector< vector<float> > Pintail::calcExpected(int start, int end) {
	try {
		
		vector< vector<float> > temp; temp.resize(end-start);
		
		//for each sequence
		int count = 0;
		for(int i = start; i < end; i++){
			
			float coef = seqCoef[i];
			
			//for each window
			vector<float> queryExpected;
			for (int m = 0; m < windows[i].size(); m++) {		
				float expected = Qav[i][m] * coef;
				queryExpected.push_back(expected);	
//cout << "average variabilty over window = " << averageProbability[m] << " coef = " << coef << " ei = "  << expected << '\t' <<  "window = " << m << endl;
			}
			
			temp[count] = queryExpected;
			
			count++;
			
		}
		
		return temp;
				
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "calcExpected");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float> Pintail::calcDE(int start, int end) {
	try {
		
		vector<float> temp; temp.resize(end-start);
	
		//for each sequence
		int count = 0;
		for(int i = start; i < end; i++){
			
			vector<float> obs = obsDistance[i];
			vector<float> exp = expectedDistance[i];
			
//	cout << "difference between obs and exp = " << abs(obs[m] - exp[m]) << endl;	
			//for each window
			float sum = 0.0;  //sum = sum from 1 to m of (oi-ei)^2
			for (int m = 0; m < windows[i].size(); m++) { 		sum += ((obs[m] - exp[m]) * (obs[m] - exp[m]));		}
			
			float de = sqrt((sum / (windows[i].size() - 1)));
			
			temp[count] = de;
			count++;
		}
		
		return temp;
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
			//if ( (seqs.size() - gaps) == 0 ) {  highFreq = 1.0;  }			
			//else { highFreq = highest / (float) (seqs.size() - gaps);	 }
			highFreq = highest / (float) seqs.size();
cout << i << '\t' << highFreq << endl;
			
			float Pi;
			Pi =  (highFreq - 0.25) / 0.75; 
			
			//saves this for later
			outFreq << i << '\t' << Pi << endl;
				
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
vector< vector<float> > Pintail::findQav(int start, int end) {
	try {
		vector< vector<float> > averages; 
		map<int, int>::iterator it;
		
		for(int i = start; i < end; i++){
		
			//for each window find average
			vector<float> temp;
			for (int m = 0; m < windows[i].size(); m++) {
				
				float average = 0.0;
				
				it = trim[i].begin();  //trim[i] is a map of where this sequence was trimmed
				
				//while you are in the window for this sequence
				for (int j = windows[i][m]+it->first; j < (windows[i][m]+windowSizes[i]); j++) {   average += probabilityProfile[j];	}
				
				average = average / windowSizes[i];
	//cout << average << endl;			
				//save this windows average
				temp.push_back(average);
			}
			
			//save this qav
			averages.push_back(temp);
		}
		
		return averages;
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "findQav");
		exit(1);
	}
}
//***************************************************************************************************************
vector<float> Pintail::getCoef(int start, int end) {
	try {
		vector<float> coefs;
		coefs.resize(end-start);
		
		//find a coef for each sequence
		int count = 0;
		for(int i = start; i < end; i++){
		
			//find average prob for this seqs windows
			float probAverage = 0.0;
			for (int j = 0; j < Qav[i].size(); j++) {   probAverage += Qav[i][j];	}
			probAverage = probAverage / (float) Qav[i].size();
	cout << "(sum of ai) / m = " << probAverage << endl;		

			vector<float> temp = obsDistance[i];
			
			//find observed average 
			float obsAverage = 0.0;
			for (int j = 0; j < temp.size(); j++) {   obsAverage += temp[j];	}
			obsAverage = obsAverage / (float) temp.size();
cout << "(sum of oi) / m = " << obsAverage << endl;		
			float coef = obsAverage / probAverage;
		
			//save this sequences coefficient
			coefs[count] = coef;
			
			count++;
		}
		
						
		return coefs;
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
		vector< map <int, int> > t; t.resize(querySeqs.size());
		
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
					t[i] = trim[i];
					
					count++;
				}
				
				
				
				vector< vector<int> > temp = findWindows(lines[process]->start, lines[process]->end);
				
				//move into best
				count = 0;
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					win[i] = temp[count];
					count++;
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
		trim = t;
#else
		windows = findWindows(lines[0]->start, lines[0]->end);

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
		
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
				
				vector< vector<float> > temp;
				vector<float> tempde;
				int count = 0;
				
				
				temp = calcObserved(lines[process]->start, lines[process]->end);
				count = 0;
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					obs[i] = temp[count];
					count++;
				}

				temp = findQav(lines[process]->start, lines[process]->end);
				count = 0;
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					Qav[i] = temp[count];
					count++;
				}
				
				tempde = getCoef(lines[process]->start, lines[process]->end);
				count = 0;
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					seqCoef[i] = tempde[count];
					count++;
				}
				
				temp = calcExpected(lines[process]->start, lines[process]->end);
				count = 0;
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					exp[i] = temp[count];
					count++;
				}

				
				tempde = calcDE(lines[process]->start, lines[process]->end); 
				count = 0;
				for (int i = lines[process]->start; i < lines[process]->end; i++) {
					de[i] = tempde[count];
					count++;
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
		
#else
		bestfit = findPairs(lines[0]->start, lines[0]->end);
		obsDistance = calcObserved(lines[0]->start, lines[0]->end);
		Qav = findQav(lines[0]->start, lines[0]->end);
		seqCoef = getCoef(lines[0]->start, lines[0]->end);
		expectedDistance = calcExpected(lines[0]->start, lines[0]->end);
		DE = calcDE(lines[0]->start, lines[0]->end); 

#endif		
	}
	catch(exception& e) {
		errorOut(e, "Pintail", "createProcesses");
		exit(1);
	}
}

//***************************************************************************************************************


