/*
 *  alignedsimilarity.cpp
 *  Mothur
 *
 *  Created by westcott on 8/18/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "alignedsimilarity.h"
#include "ignoregaps.h"


//***************************************************************************************************************
AlignSim::AlignSim(string filename, string temp) {  fastafile = filename;  templateFile = temp;  }
//***************************************************************************************************************

AlignSim::~AlignSim() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)		{  delete querySeqs[i];		}
		for (int i = 0; i < templateSeqs.size(); i++)	{  delete templateSeqs[i];	}
		delete distCalc;
	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "~AlignSim");
		exit(1);
	}
}	
//***************************************************************************************************************
void AlignSim::print(ostream& out) {
	try {
		
		mothurOutEndLine();
		
		for (int i = 0; i < querySeqs.size(); i++) {
			
			int j = 0;  float largest = -10;
			//find largest sim value
			for (int k = 0; k < IS[i].size(); k++) {
				//is this score larger
				if (IS[i][k].score > largest) {
					j = k;
					largest = IS[i][k].score;
				}
			}
			
			//find parental similarity
			distCalc->calcDist(*(IS[i][j].leftParent), *(IS[i][j].rightParent));
			float dist = distCalc->getDist();
			
			//convert to similarity
			dist = (1 - dist) * 100;

			//warn about parental similarity - if its above 82% may not detect a chimera 
			if (dist >= 82) { mothurOut("When the chimeras parental similarity is above 82%, detection rates drop signifigantly.");  mothurOutEndLine(); }
			
			int index = ceil(dist);
			
			if (index == 0) { index=1;  }
	
			//is your DE value higher than the 95%
			string chimera;
			if (IS[i][j].score > quantile[index-1][4])		{	chimera = "Yes";	}
			else										{	chimera = "No";		}			
			
			out << querySeqs[i]->getName() <<  "\tparental similarity: " << dist << "\tIS: " << IS[i][j].score << "\tbreakpoint: " << IS[i][j].midpoint << "\tchimera flag: " << chimera << endl;
			
			if (chimera == "Yes") {
				mothurOut(querySeqs[i]->getName() + "\tparental similarity: " + toString(dist) + "\tIS: " + toString(IS[i][j].score) + "\tbreakpoint: " + toString(IS[i][j].midpoint) + "\tchimera flag: " + chimera); mothurOutEndLine();
			}
			out << "Improvement Score\t";
			
			for (int r = 0; r < IS[i].size(); r++) {  out << IS[i][r].score << '\t';  }
			out << endl;
		}
	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "print");
		exit(1);
	}
}

//***************************************************************************************************************
void AlignSim::getChimeras() {
	try {
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		templateSeqs = readSeqs(templateFile);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		IS.resize(numSeqs);
		
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
		
		//find window breaks
		windowBreak = findWindows();
//for (int i = 0; i < windowBreak.size(); i++) { cout << windowBreak[i] << '\t';  }
//cout << endl;

		mothurOut("Finding the IS values for your sequences..."); cout.flush();
		if (processors == 1) { 
			IS = findIS(lines[0]->start, lines[0]->end, querySeqs);
		}else {		IS = createProcessesIS(querySeqs, lines);		}
		mothurOut("Done."); mothurOutEndLine();
		
		//quantiles are used to determine whether the de values found indicate a chimera
		if (quanfile != "") {  quantile = readQuantiles();  }
		else {
			
			mothurOut("Calculating quantiles for your template.  This can take a while...  I will output the quantiles to a .alignedsimilarity.quan file that you can input them using the quantiles parameter next time you run this command.  Providing the .quan file will dramatically improve speed.    "); cout.flush();
			//get IS quantiles as a reference to these IS scores
			//you are assuming the template is free of chimeras and therefore will give you a baseline as to the scores you would like to see
			if (processors == 1) { 
				templateIS = findIS(templateLines[0]->start, templateLines[0]->end, templateSeqs);
			}else {		templateIS = createProcessesIS(templateSeqs, templateLines);		}
//cout << "here" << endl;			
			ofstream out4;
			string o;
			
			o = getRootName(templateFile) + "alignedsimilarity.quan";
			
			openOutputFile(o, out4);
			
			//adjust quantiles
			//construct table
			quantile.resize(100);
		
			for (int i = 0; i < templateIS.size(); i++) {
				
				if (templateIS[i].size() != 0) {
					int j = 0; float score = -1000;
					//find highest IS score
	//cout << templateIS[i].size() << endl;
					for (int k = 0; k < templateIS[i].size(); k++) {
						if (templateIS[i][k].score > score) {
							score = templateIS[i][k].score;
							j = k;
						}
					}
	//cout << j << endl;		
					//find similarity of parents
					distCalc->calcDist(*(templateIS[i][j].leftParent), *(templateIS[i][j].rightParent));
					float dist = distCalc->getDist();
			
					//convert to similarity
					dist = (1 - dist) * 100;
	//	cout << dist << endl;	
					int index = ceil(dist);
	//	cout << "index = " << index << endl;	
					if (index == 0) {	index = 1;	}
					quantile[index-1].push_back(templateIS[i][j].score);
				}
			}
		
		
			for (int i = 0; i < quantile.size(); i++) {
				vector<float> temp;
				
				if (quantile[i].size() == 0) {
					//in case this is not a distance found in your template files
					for (int g = 0; g < 6; g++) {
						temp.push_back(0.0);
					}
				}else{
					
					sort(quantile[i].begin(), quantile[i].end());
					
					//save 10%
					temp.push_back(quantile[i][int(quantile[i].size() * 0.10)]);
					//save 25%
					temp.push_back(quantile[i][int(quantile[i].size() * 0.25)]);
					//save 50%
					temp.push_back(quantile[i][int(quantile[i].size() * 0.5)]);
					//save 75%
					temp.push_back(quantile[i][int(quantile[i].size() * 0.75)]);
					//save 95%
					temp.push_back(quantile[i][int(quantile[i].size() * 0.95)]);
					//save 99%
					temp.push_back(quantile[i][int(quantile[i].size() * 0.99)]);
					
				}
				
				//output quan value
				out4 << i+1 << '\t';				
				for (int u = 0; u < temp.size(); u++) {   out4 << temp[u] << '\t'; }
				out4 << endl;
				
				quantile[i] = temp;
				
			}
			
			out4.close();
			
			mothurOut("Done."); mothurOutEndLine();

		}
		
		//free memory
		for (int i = 0; i < lines.size(); i++)					{	delete lines[i];				}
		for (int i = 0; i < templateLines.size(); i++)			{	delete templateLines[i];		}
			
	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
vector<int> AlignSim::findWindows() {
	try {
		
		vector<int> win; 
		
		if (increment > querySeqs[0]->getAligned().length()) {  mothurOut("You have selected an increment larger than the length of your sequences.  I will use the default of 25.");  increment = 25; }
		
		for (int m = increment;  m < (querySeqs[0]->getAligned().length() - increment); m+=increment) {  win.push_back(m);  }

		return win;
	
	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "findWindows");
		exit(1);
	}
}

//***************************************************************************************************************
vector< vector<sim>  > AlignSim::findIS(int start, int end, vector<Sequence*> seqs) {
	try {
		
		vector< vector<sim> >  isValues;
		isValues.resize(seqs.size());
		
		//for each sequence
		for(int i = start; i < end; i++){
			
			vector<sim> temp;  temp.resize(windowBreak.size());
			
			mothurOut("Finding IS value for sequence " + toString(i)); mothurOutEndLine();
			
			//for each window
			for (int m = 0; m < windowBreak.size(); m++) {
				
				vector<Sequence*> closest;  //left, right, overall
				vector<int> similarity; //left+right and overall
				
				//find closest left, closest right and closest overall
				closest = findClosestSides(seqs[i], windowBreak[m], similarity, i);
				
				int totalNumBases = seqs[i]->getUnaligned().length();
				
				//IS = left+right-overall
				float is = ((similarity[0]+similarity[1]) / (float) totalNumBases) - (similarity[2] / (float) totalNumBases);
				
				///save IS, leftparent, rightparent, breakpoint
				temp[m].leftParent = closest[0];
				temp[m].rightParent = closest[1];
				temp[m].score = is;
				temp[m].midpoint = windowBreak[m];
//cout << is << '\t';
			}
//cout << endl;			
			isValues[i] = temp;
		
		}
		
		return isValues;
	
	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "findIS");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence*> AlignSim::findClosestSides(Sequence* seq, int breakpoint, vector<int>& sim, int i) {
	try{
	
		vector<Sequence*> closest;
		
		Sequence query, queryLeft, queryRight;
		string frag = seq->getAligned();
		string fragLeft = frag.substr(0, breakpoint);
		string fragRight = frag.substr(breakpoint, frag.length());
		
		//get pieces
		query = *(seq);  
		queryLeft = *(seq);  
		queryRight = *(seq);
		queryLeft.setAligned(fragLeft);
		queryRight.setAligned(fragRight);
		
		//initialize
		Sequence* overall = templateSeqs[0];
		Sequence* left = templateSeqs[0];
		Sequence* right = templateSeqs[0];
		
		float smallestOverall, smallestLeft, smallestRight;
		smallestOverall = 1000;  smallestLeft = 1000;  smallestRight = 1000;
		
		//go through the templateSeqs and search for the closest
		for(int j = 0; j < templateSeqs.size(); j++){
			
			//so you don't pick yourself
			if (j != i) {
				Sequence temp, tempLeft, tempRight;
				string fragTemp = templateSeqs[j]->getAligned();
				string fragTempLeft = fragTemp.substr(0, breakpoint);
				string fragTempRight = fragTemp.substr(breakpoint, fragTemp.length());
				
				//get pieces
				temp = *(templateSeqs[j]); 
				tempLeft = *(templateSeqs[j]); 
				tempRight = *(templateSeqs[j]);
				tempLeft.setAligned(fragTempLeft);
				tempRight.setAligned(fragTempRight);
				
				//find overall dist
				distCalc->calcDist(query, temp);
				float dist = distCalc->getDist();	
				
				if (dist < smallestOverall) { 
					overall = templateSeqs[j];
					smallestOverall = dist;
				}
				
				//find left dist
				distCalc->calcDist(queryLeft, tempLeft);
				dist = distCalc->getDist();
				
				if (dist < smallestLeft) { 
					left = templateSeqs[j];
					smallestLeft = dist;
				}
				
				//find left dist
				distCalc->calcDist(queryRight, tempRight);
				dist = distCalc->getDist();
				
				if (dist < smallestRight) { 
					right = templateSeqs[j];
					smallestRight = dist;
				}
			}

		}
			
		closest.push_back(left);
		closest.push_back(right);
		closest.push_back(overall);
		
		//fill sim with number of matched bases
		Sequence tempL = *(left);
		string tempLFrag = tempL.getAligned();
		tempL.setAligned(tempLFrag.substr(0, breakpoint));
		
		int bothMatches = findNumMatchedBases(queryLeft, tempL);
		sim.push_back(bothMatches);

		Sequence tempR = *(right);
		string tempRFrag = tempR.getAligned();
		tempR.setAligned(tempRFrag.substr(breakpoint, tempRFrag.length()));
		
		bothMatches = findNumMatchedBases(queryRight, tempR);
		sim.push_back(bothMatches);

		bothMatches = findNumMatchedBases(query, *(overall));
		sim.push_back(bothMatches);
			
		return closest;


	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "findClosestSides");
		exit(1);
	}
}
//***************************************************************************************************************

int AlignSim::findNumMatchedBases(Sequence seq, Sequence temp) {
	try{
	
		int num = 0;
		
		string query = seq.getAligned();
		string subject = temp.getAligned();
//cout << seq.getName() << endl << endl;
//cout << temp.getName()  << endl << endl;
		
		for (int i = 0; i < query.length(); i++) {
			//count matches if query[i] is a base and subject is same bases or gap
			if (isalpha(query[i])) {
				if(!isalpha(subject[i])) {  num++;  } 
				else if (query[i] == subject[i])  { num++;  }
				else {  }
			}
		}
//cout << "num = " << num << endl;		
		return num;
	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "findNumMatchedBases");
		exit(1);
	}
}

/**************************************************************************************************/
vector< vector<sim> > AlignSim::createProcessesIS(vector<Sequence*> seqs, vector<linePair*> linesToProcess) {
	try {
	vector< vector<sim> > localIs; localIs.resize(seqs.size());
	
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
				
				mothurOut("Finding IS values for sequences " + toString(linesToProcess[process]->start) + " to " + toString(linesToProcess[process]->end)); mothurOutEndLine();
				localIs = findIS(linesToProcess[process]->start, linesToProcess[process]->end, seqs);
				mothurOut("Done finding IS values for sequences " +  toString(linesToProcess[process]->start) + " to " + toString(linesToProcess[process]->end)); mothurOutEndLine();
				
				//write out data to file so parent can read it
				ofstream out;
				string s = toString(getpid()) + ".temp";
				openOutputFile(s, out);
				
				//output pairs
				for (int i = linesToProcess[process]->start; i < linesToProcess[process]->end; i++) {
					 out << localIs[i].size() << endl;
					 for (int j = 0; j < localIs[i].size(); j++) {
						localIs[i][j].leftParent->printSequence(out);
						localIs[i][j].rightParent->printSequence(out);
						out << ">" << '\t' << localIs[i][j].score << '\t' << localIs[i][j].midpoint << endl;
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
			for (int k = linesToProcess[i]->start; k < linesToProcess[i]->end; k++) {
				int size;
				in >> size;
				gobble(in);
				
				vector<sim> tempVector;
				
				for (int j = 0; j < size; j++) {
				
					sim temp;
					
					temp.leftParent = new Sequence(in);
					gobble(in);
	//temp.leftParent->printSequence(cout);				
					temp.rightParent = new Sequence(in);
					gobble(in);
//temp.rightParent->printSequence(cout);
					string throwaway;
					in >> throwaway >> temp.score >> temp.midpoint;
	//cout << temp.score << '\t' << temp.midpoint << endl;				
					gobble(in);
					
					tempVector.push_back(temp);
				}
				
				localIs[k] = tempVector;
			}
			
			in.close();
			remove(s.c_str());
		}
			
	
#else
		localIs = findIS(linesToProcess[0]->start, linesToProcess[0]->end, seqs);
#endif	

		return localIs;
	}
	catch(exception& e) {
		errorOut(e, "AlignSim", "createProcessesIS");
		exit(1);
	}
}

//***************************************************************************************************************
