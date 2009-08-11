/*
 *  mallard.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 8/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "mallard.h"

//***************************************************************************************************************

Mallard::Mallard(string filename) {  fastafile = filename;  }
//***************************************************************************************************************

Mallard::~Mallard() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)		{  delete querySeqs[i];		}
	}
	catch(exception& e) {
		errorOut(e, "Mallard", "~Mallard");
		exit(1);
	}
}	
//***************************************************************************************************************
void Mallard::print(ostream& out) {
	try {
		
		for (int i = 0; i < querySeqs.size(); i++) {
			
			out << querySeqs[i]->getName() << "\thighest de value = " << highestDE[i] << "\tpenalty score = " << marked[i] << endl;
			cout << querySeqs[i]->getName() << "\tpenalty score = " << marked[i] << endl;
						
		}
	}
	catch(exception& e) {
		errorOut(e, "Mallard", "print");
		exit(1);
	}
}

//***************************************************************************************************************
void Mallard::getChimeras() {
	try {
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		windowSizes.resize(numSeqs, window);
		quantilesMembers.resize(100);  //one for every percent mismatch
		highestDE.resize(numSeqs, 0.0);  //contains the highest de value for each seq
		
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
			
		#else
			lines.push_back(new linePair(0, numSeqs));
		#endif
		
		decalc = new DeCalculator();
		
		//if the user does enter a mask then you want to keep all the spots in the alignment
		if (seqMask.length() == 0)	{	decalc->setAlignmentLength(querySeqs[0]->getAligned().length());	}
		else						{	decalc->setAlignmentLength(seqMask.length());						}
		
		decalc->setMask(seqMask);
				
		//find P
		mothurOut("Getting conservation... "); cout.flush();
		probabilityProfile = decalc->calcFreq(querySeqs, fastafile); 
	
		//make P into Q
		for (int i = 0; i < probabilityProfile.size(); i++)  {	probabilityProfile[i] = 1 - probabilityProfile[i];  } 
		mothurOut("Done."); mothurOutEndLine();
		
		//mask sequences if the user wants to 
		if (seqMask != "") {
			//mask querys
			for (int i = 0; i < querySeqs.size(); i++) {
				decalc->runMask(querySeqs[i]);
			}
		}
						
		mothurOut("Calculating DE values..."); cout.flush();
		if (processors == 1) { 
			quantilesMembers = decalc->getQuantiles(querySeqs, windowSizes, window, probabilityProfile, increment, 0, querySeqs.size(), highestDE);
		}else {		createProcessesQuan();		}
		mothurOut("Done."); mothurOutEndLine();
		
		mothurOut("Ranking outliers..."); cout.flush();
		marked = decalc->returnObviousOutliers(quantilesMembers, querySeqs.size());
		mothurOut("Done."); mothurOutEndLine();

	
		//free memory
		for (int i = 0; i < lines.size(); i++)					{	delete lines[i];		}		
		delete decalc;
	}
	catch(exception& e) {
		errorOut(e, "Mallard", "getChimeras");
		exit(1);
	}
}
/**************************************************************************************************/

void Mallard::createProcessesQuan() {
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
				
				quantilesMembers = decalc->getQuantiles(querySeqs, windowSizes, window, probabilityProfile, increment, lines[process]->start, lines[process]->end, highestDE);
				
				//write out data to file so parent can read it
				ofstream out;
				string s = toString(getpid()) + ".temp";
				openOutputFile(s, out);
				
								
				//output observed distances
				for (int i = 0; i < quantilesMembers.size(); i++) {
					out << quantilesMembers[i].size() << '\t';
					for (int j = 0; j < quantilesMembers[i].size(); j++) {
						out << quantilesMembers[i][j].score << '\t' << quantilesMembers[i][j].member1 << '\t' << quantilesMembers[i][j].member2 << '\t';
					}
					out << endl;
				}
				
				out << highestDE.size() << endl;
				for (int i = 0; i < highestDE.size(); i++) {
					out << highestDE[i] << '\t';
				}
				out << endl;
				
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
			
			vector< vector<quanMember> > quan; 
			quan.resize(100);
			
			//get quantiles
			for (int m = 0; m < quan.size(); m++) {
				int num;
				in >> num; 
				
				gobble(in);

				vector<quanMember> q;  float w; int b, n;
				for (int j = 0; j < num; j++) {
					in >> w >> b >> n;
	//cout << w << '\t' << b << '\t' n << endl;
					quanMember newMember(w, b, n);
					q.push_back(newMember);
				}
//cout << "here" << endl;
				quan[m] = q;
//cout << "now here" << endl;
				gobble(in);
			}
			
	
			//save quan in quantiles
			for (int j = 0; j < quan.size(); j++) {
				//put all values of q[i] into quan[i]
				for (int l = 0; l < quan[j].size(); l++) {  quantilesMembers[j].push_back(quan[j][l]);   }
				//quantilesMembers[j].insert(quantilesMembers[j].begin(), quan[j].begin(), quan[j].end());
			}
			
			int num;
			in >> num;  gobble(in);
			
			int count = lines[process]->start;
			for (int s = 0; s < num; s++) {
				float high;
				in >> high;
				
				highestDE[count] = high;
				count++;
			}
			
			in.close();
			remove(s.c_str());
		}
		
#else
		quantilesMembers = decalc->getQuantiles(querySeqs, windowSizes, window, probabilityProfile, increment, 0, querySeqs.size(), highestDE);
#endif		
	}
	catch(exception& e) {
		errorOut(e, "Mallard", "createProcessesQuan");
		exit(1);
	}
}
//***************************************************************************************************************


