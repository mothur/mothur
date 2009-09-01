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
	
		distCalc = new eachGapDist();
		decalc = new DeCalculator();
		
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

		//mask sequences if the user wants to 
		if (seqMask != "") {
			//mask querys
			for (int i = 0; i < querySeqs.size(); i++) {
				decalc->runMask(querySeqs[i]);
			}
		
			//mask templates
			for (int i = 0; i < templateSeqs.size(); i++) {
				decalc->runMask(templateSeqs[i]);
			}
		}
		
		if (filter) {
			vector<Sequence*> temp = templateSeqs;
			for (int i = 0; i < querySeqs.size(); i++) { temp.push_back(querySeqs[i]);  }
			
			createFilter(temp);
			
			runFilter(querySeqs);
			runFilter(templateSeqs);
		}

		//trim sequences - this follows ccodes remove_extra_gaps 
		//just need to pass it query and template since closest points to template
		trimSequences();
		
		//windows are equivalent to words - ccode paper recommends windows are between 5% and 20% on alignment length().  
		//Our default will be 10% and we will warn if user tries to use a window above or below these recommendations
		windows = findWindows();  
		
		//remove sequences that are more than 20% different and less than 0.5% different - may want to allow user to specify this later
		for (int i = 0; i < closest.size(); i++) {
			removeBadReferenceSeqs(closest[i], i);
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
/***************************************************************************************************************/
//ccode algo says it does this to "Removes the initial and final gaps to avoid biases due to incomplete sequences."
void Ccode::trimSequences() {
	try {
		
		int frontPos = 0;  //should contain first position in all seqs that is not a gap character
		int rearPos = querySeqs[0]->getAligned().length();
		
		//********find first position in all seqs that is a non gap character***********//
		//find first position all query seqs that is a non gap character
		for (int i = 0; i < querySeqs.size(); i++) {
			
			string aligned = querySeqs[i]->getAligned();
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
		
		//find first position all template seqs that is a non gap character
		for (int i = 0; i < templateSeqs.size(); i++) {
			
			string aligned = templateSeqs[i]->getAligned();
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

		
		//********find last position in all seqs that is a non gap character***********//
		//find last position all query seqs that is a non gap character
		for (int i = 0; i < querySeqs.size(); i++) {
			
			string aligned = querySeqs[i]->getAligned();
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
		
		//find last position all template seqs that is a non gap character
		for (int i = 0; i < templateSeqs.size(); i++) {
			
			string aligned = templateSeqs[i]->getAligned();
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

		
		//check to make sure that is not whole seq
		if ((rearPos - frontPos - 1) <= 0) {  mothurOut("Error, when I trim your sequences, the entire sequence is trimmed."); mothurOutEndLine(); exit(1);  }
		
		//***********trim all seqs to that position***************//
		for (int i = 0; i < querySeqs.size(); i++) {
			
			string aligned = querySeqs[i]->getAligned();
			
			//between the two points
			aligned = aligned.substr(frontPos, (rearPos-frontPos-1));
			
			querySeqs[i]->setAligned(aligned);
		}
		
		for (int i = 0; i < templateSeqs.size(); i++) {
			
			string aligned = templateSeqs[i]->getAligned();
			
			//between the two points
			aligned = aligned.substr(frontPos, (rearPos-frontPos-1));
			
			templateSeqs[i]->setAligned(aligned);
		}
	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "trimSequences");
		exit(1);
	}

}
/***************************************************************************************************************/
vector<int> Ccode::findWindows() {
	try {
		
		vector<int> win; 
		int length = querySeqs[0]->getAligned().length();
		
		//default is wanted = 10% of total length
		if (window > length) { 
			mothurOut("You have slected a window larger than your sequence length after all filters, masks and trims have been done. I will use the default 10% of sequence length.");
			window = length / 10;
		}else if (window == 0) { window = length / 10;  }
		else if (window > (length / 20)) {
			mothurOut("You have selected a window that is larger than 20% of your sequence length.  This is not recommended, but I will continue anyway."); mothurOutEndLine();
		}else if (window < (length / 5)) {
			mothurOut("You have selected a window that is smaller than 5% of your sequence length.  This is not recommended, but I will continue anyway."); mothurOutEndLine();
		}
		
		//save starting points of each window
		for (int m = 0;  m < (length-window); m+=window) {  win.push_back(m);  }

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
				if (seqA[i] != seqB[i]) { numDiff++; }
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
void Ccode::removeBadReferenceSeqs(vector<Sequence*>& seqs, int query) {
	try {
		
		vector< vector<int> > numDiffBases;
		numDiffBases.resize(seqs.size());
		//initialize to 0
		for (int i = 0; i < numDiffBases.size(); i++) { numDiffBases[i].resize(seqs.size(),0); }
		
		int length = seqs[0]->getAligned().length();
		
		//calc differences from each sequence to everyother seq in the set
		for (int i = 0; i < seqs.size(); i++) {
			
			string seqA = seqs[i]->getAligned();
			
			//so you don't calc i to j and j to i since they are the same
			for (int j = 0; j < i; j++) {
				
				string seqB = seqs[j]->getAligned();
				
				//compare strings
				int numDiff = getDiff(seqA, seqB);
				
				numDiffBases[i][j] = numDiff;
				numDiffBases[j][i] = numDiff;
			}
		}
		
		//initailize remove to 0
		vector<int> remove;  remove.resize(seqs.size(), 0);
		
		//check each numDiffBases and if any are higher than threshold set remove to 1 so you can remove those seqs from the closest set
		for (int i = 0; i < numDiffBases.size(); i++) {
			for (int j = 0; j < numDiffBases[i].size(); j++) {
				
				//are you more than 20% different
				if (numDiffBases[i][j] > ((20*length) / 100))		{  remove[j] = 1;  }
				//are you less than 0.5% different
				if (numDiffBases[i][j] < ((0.5*length) / 100))	{  remove[j] = 1;  }
			}
		}
		
		int numSeqsLeft = 0;
		
		//count seqs that are not going to be removed
		for (int i = 0; i < remove.size(); i++) {  
			if (remove[i] == 0)  { numSeqsLeft++;  }
		}
		
		//if you have enough then remove bad ones
		if (numSeqsLeft >= 3) {
			vector<Sequence*> goodSeqs;
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
vector<float> Ccode::getAverageRef(vector<Sequence*> ref) {
	try {
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getAverageRef");
		exit(1);
	}
}
/**************************************************************************************************/
vector<float> Ccode::getAverageQuery (vector<Sequence*> ref, int query) {
	try {
	
	
	}
	catch(exception& e) {
		errorOut(e, "Ccode", "getAverageQuery");
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
				
				vector<Sequence*> tempVector;
				
				for (int j = 0; j < numWanted; j++) {
				
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

