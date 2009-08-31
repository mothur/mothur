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


//***************************************************************************************************************
Ccode::Ccode(string filename, string temp) {  fastafile = filename;  templateFile = temp;  }
//***************************************************************************************************************

Ccode::~Ccode() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)		{  delete querySeqs[i];		}
		for (int i = 0; i < templateSeqs.size(); i++)	{  delete templateSeqs[i];	}
		delete distCalc;
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
		
		for (int i = 0; i < querySeqs.size(); i++) {
			
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
			
			//find breakup of templatefile for quantiles
			if (processors == 1) {   templateLines.push_back(new linePair(0, templateSeqs.size()));  }
			else { 
				for (int i = 0; i < processors; i++) {
					templateLines.push_back(new linePair());
					templateLines[i]->start = int (sqrt(float(i)/float(processors)) * templateSeqs.size());
					templateLines[i]->end = int (sqrt(float(i+1)/float(processors)) * templateSeqs.size());
				}
			}
		#else
			lines.push_back(new linePair(0, numSeqs));
			templateLines.push_back(new linePair(0, templateSeqs.size()));
		#endif
	
		distCalc = new ignoreGaps();
		
		//find closest
		if (processors == 1) { 
			mothurOut("Finding top matches for sequences... "); cout.flush();
			closest = findClosest(lines[0]->start, lines[0]->end, numWanted);
			mothurOut("Done."); mothurOutEndLine();
		}else {		createProcessesClosest();		}

		
		for (int i = 0; i < closest.size(); i++) {
			cout << querySeqs[i]->getName() << ": ";
			for (int j = 0; j < closest[i].size(); j++) {
			
				cout << closest[i][j]->getName() << '\t';
			}
			cout << endl;
		}
			
		//free memory
		for (int i = 0; i < lines.size(); i++)					{	delete lines[i];				}
		for (int i = 0; i < templateLines.size(); i++)			{	delete templateLines[i];		}
			
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getChimeras");
		exit(1);
	}
}
/***************************************************************************************************************
vector<int> Ccode::findWindows() {
	try {
		
		vector<int> win; 
		
		if (increment > querySeqs[0]->getAligned().length()) {  mothurOut("You have selected an increment larger than the length of your sequences.  I will use the default of 25.");  increment = 25; }
		
		for (int m = increment;  m < (querySeqs[0]->getAligned().length() - increment); m+=increment) {  win.push_back(m);  }

		return win;
	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "findWindows");
		exit(1);
	}
}
*/
//***************************************************************************************************************
vector< vector<Sequence*> > Ccode::findClosest(int start, int end, int numWanted) {
	try{
	
		vector< vector<Sequence*> > topMatches;  topMatches.resize(querySeqs.size());
	
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
				topMatches[j].push_back(distances[h].seq);
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
					 out << closest[i].size() << endl;
					 for (int j = 0; j < closest[i].size(); j++) {
						closest[i][j]->printSequence(out);
					 }
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
				int size;
				in >> size;
				gobble(in);
				
				vector<Sequence*> tempVector;
				
				for (int j = 0; j < size; j++) {
				
					Sequence* temp = new Sequence(in);
					gobble(in);
						
					tempVector.push_back(temp);
				}
				
				closest[k] = tempVector;
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

