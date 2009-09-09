/*
 *  chimeracheckrdp.cpp
 *  Mothur
 *
 *  Created by westcott on 9/8/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "chimeracheckrdp.h"
		
//***************************************************************************************************************
ChimeraCheckRDP::ChimeraCheckRDP(string filename, string temp) {  fastafile = filename;  templateFile = temp;  }
//***************************************************************************************************************

ChimeraCheckRDP::~ChimeraCheckRDP() {
	try {
		for (int i = 0; i < querySeqs.size(); i++)		{  delete querySeqs[i];		}
		delete templateDB;
		delete kmer;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraCheckRDP", "~AlignSim");
		exit(1);
	}
}	
//***************************************************************************************************************
void ChimeraCheckRDP::print(ostream& out) {
	try {
		
		mothurOutEndLine();
	/*	
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
		}*/
	}
	catch(exception& e) {
		errorOut(e, "ChimeraCheckRDP", "print");
		exit(1);
	}
}

//***************************************************************************************************************
void ChimeraCheckRDP::getChimeras() {
	try {
		
		//read in query sequences and subject sequences
		mothurOut("Reading sequences and template file... "); cout.flush();
		querySeqs = readSeqs(fastafile);
		//templateSeqs = readSeqs(templateFile);
		templateDB = new KmerDB(templateFile, kmerSize);
		mothurOut("Done."); mothurOutEndLine();
		
		int numSeqs = querySeqs.size();
		
		IS.resize(numSeqs);
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
			
		#else
			lines.push_back(new linePair(0, numSeqs));
		#endif
		
		kmer = new Kmer(kmerSize);
		
		//find closest seq to each querySeq
		for (int i = 0; i < querySeqs.size(); i++) {
			closest[i] = templateDB->findClosestSequence(querySeqs[i]);  
		}
		
		//fill seqKmerInfo for query seqs
		for (int i = 0; i < querySeqs.size(); i++) {
			seqKmerInfo[querySeqs[i]->getName()] = kmer->getKmerCounts(querySeqs[i]->getUnaligned());
		}
		
		//fill seqKmerInfo for closest
		for (int i = 0; i < closest.size(); i++) {
			seqKmerInfo[closest[i].getName()] = kmer->getKmerCounts(closest[i].getUnaligned());
		}

		
		//for each query find IS value - this should be paralellized, 
		//but paralellizing may cause you to have to recalculate some seqKmerInfo since the separate processes don't share memory after they split
		for (int i = 0; i < querySeqs.size(); i++) {
			IS[i] = findIS(i);  //fills seqKmerInfo
		}
		
		
		//free memory
		for (int i = 0; i < lines.size(); i++)					{	delete lines[i];	}
	
			
	}
	catch(exception& e) {
		errorOut(e, "ChimeraCheckRDP", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
vector<sim> ChimeraCheckRDP::findIS(int query) {
	try {
		
		vector<sim>  isValues;
		string queryName = querySeqs[query]->getName();
		string seq = querySeqs[query]->getUnaligned();
		
		mothurOut("Finding IS values for sequence " + query); mothurOutEndLine();
		
		//find total kmers you have in common with closest[query] by looking at the last entry in the vector of maps for each
		int nTotal = calcKmers(seqKmerInfo[queryName][(seqKmerInfo[queryName].size()-1)], seqKmerInfo[closest[query].getName()][(seqKmerInfo[closest[query].getName()].size()-1)]);
		
		//you don't want the starting point to be virtually at hte end so move it in 10%
		int start = seq.length() / 10;
			
		//for each window
		for (int m = start; m < (seq.length() - start); m+=increment) {
			
			sim temp;
			
			string fragLeft = seq.substr(0, m);  //left side of breakpoint
			string fragRight = seq.substr(m, seq.length());  //right side of breakpoint
			
			//make a sequence of the left side and right side
			Sequence* left = new Sequence(queryName, fragLeft);
			Sequence* right = new Sequence(queryName, fragRight);
			
			//find seqs closest to each fragment
			Sequence closestLeft = templateDB->findClosestSequence(left); 
			Sequence closestRight = templateDB->findClosestSequence(right); 
			
			map<int, int>::iterator itleft;
			map<int, int>::iterator itleftclose;
			
			//get kmer in the closest seqs
			//if it's not found calc kmer info and save, otherwise use already calculated data
			//left
			it = seqKmerInfo.find(closestLeft.getName());
			if (it == seqKmerInfo.end()) {  //you have to calc it
				seqKmerInfo[closestLeft.getName()] = kmer->getKmerCounts(closestLeft.getUnaligned());
			}
			
			//right
			it = seqKmerInfo.find(closestRight.getName());
			if (it == seqKmerInfo.end()) {  //you have to calc it
				seqKmerInfo[closestRight.getName()] = kmer->getKmerCounts(closestRight.getUnaligned());
			}
			
			//right side is tricky - since the counts grow on eachother to find the correct counts of only the right side you must subtract the counts of the left side
			//iterate through left sides map to subtract the number of times you saw things before you got the the right side
			map<int, int> rightside;
			for (itleft = seqKmerInfo[queryName][m-kmerSize].begin(); itleft != seqKmerInfo[queryName][m-kmerSize].end(); itleft++) {
				int howManyTotal = seqKmerInfo[queryName][seqKmerInfo[queryName].size()-1][itleft->first];   //times that kmer was seen in total
				
				//itleft->second is times it was seen in left side, so howmanytotal - leftside should give you right side
				int howmanyright = howManyTotal - itleft->second;
				
				//if any were seen just on the right add that ammount to map
				if (howmanyright > 0) {
					rightside[itleft->first] = howmanyright;
				}
			}
			
			//iterate through left side of the seq closest to the right fragment of query to subtract the number you saw before you reached the right side of the closest right
			//this way you can get the map for just the fragment you want to compare and not hte whole sequence
			map<int, int> rightsideclose;
			for (itleftclose = seqKmerInfo[closestRight.getName()][m-kmerSize].begin(); itleftclose != seqKmerInfo[closestRight.getName()][m-kmerSize].end(); itleftclose++) {
				int howManyTotal = seqKmerInfo[closestRight.getName()][seqKmerInfo[closestRight.getName()].size()-1][itleftclose->first];   //times that kmer was seen in total
				
				//itleft->second is times it was seen in left side, so howmanytotal - leftside should give you right side
				int howmanyright = howManyTotal - itleftclose->second;
				
				//if any were seen just on the right add that ammount to map
				if (howmanyright > 0) {
					rightsideclose[itleftclose->first] = howmanyright;
				}
			}

			int nLeft = calcKmers(seqKmerInfo[closestLeft.getName()][m-kmerSize], seqKmerInfo[queryName][m-kmerSize]);
			int nRight = calcKmers(rightsideclose, rightside);
			
			int is = nLeft + nRight - nTotal;
			
			//save IS, leftparent, rightparent, breakpoint
			temp.leftParent = closestLeft.getName();
			temp.rightParent = closestRight.getName();
			temp.score = is;
			temp.midpoint = m;
			
			isValues.push_back(temp);
			
			delete left;
			delete right;
		}	
		
		return isValues;
	
	}
	catch(exception& e) {
		errorOut(e, "ChimeraCheckRDP", "findIS");
		exit(1);
	}
}

//***************************************************************************************************************
//find the smaller map and iterate through it and count kmers in common
int ChimeraCheckRDP::calcKmers(map<int, int> query, map<int, int> subject) {
	try{
		
		int common = 0;
		map<int, int>::iterator small;
		map<int, int>::iterator large;
		
		if (query.size() < subject.size()) {
		
			for (small = query.begin(); small != query.end(); small++) {
				large = subject.find(small->first);
				
				//if you found it they have that kmer in common
				if (large != subject.end()) {	common++;	}
			}
			
		}else { 
		 
			for (small = subject.begin(); small != subject.end(); small++) {
				large = query.find(small->first);
				
				//if you found it they have that kmer in common
				if (large != query.end()) {		common++;	}
			}
		}
		
		return common;
		
	}
	catch(exception& e) {
		errorOut(e, "ChimeraCheckRDP", "calcKmers");
		exit(1);
	}
}

//***************************************************************************************************************

