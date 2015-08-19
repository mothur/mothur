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
ChimeraCheckRDP::ChimeraCheckRDP(string filename, string temp, string n, bool s, int inc, int k, string o) : Chimera() { 
	try {
		fastafile = filename; 
		templateFileName = temp;  
		name = n;
		svg = s;
		increment = inc;
		kmerSize = k;
		outputDir = o; 
		
		templateDB = new AlignmentDB(templateFileName, "kmer", kmerSize, 0.0,0.0,0.0,0.0, rand());
		m->mothurOutEndLine();
		
		kmer = new Kmer(kmerSize);
		
		if (name != "") { 
			readName(name);  //fills name map with names of seqs the user wants to have .svg for.  
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "ChimeraCheckRDP");
		exit(1);
	}
}
//***************************************************************************************************************

ChimeraCheckRDP::~ChimeraCheckRDP() {
	try {
		delete templateDB;
		delete kmer;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "~ChimeraCheckRDP");
		exit(1);
	}
}	
//***************************************************************************************************************
Sequence ChimeraCheckRDP::print(ostream& out, ostream& outAcc) {
	try {
		
		m->mothurOut("Processing: " + querySeq->getName()); m->mothurOutEndLine();
		
		out << querySeq->getName() << endl;
		out << "IS scores: " << '\t';
			
		for (int k = 0; k < IS.size(); k++) {
			out << IS[k].score << '\t'; 
		}
		out << endl;
		
		if (svg) {
			if (name != "") { //if user has specific names
				map<string, string>::iterator it = names.find(querySeq->getName());
				
				if (it != names.end()) { //user wants pic of this
					makeSVGpic(IS);  //zeros out negative results
				}
			}else{//output them all
				makeSVGpic(IS);  //zeros out negative results
			}
		}
		
		return *querySeq;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "print");
		exit(1);
	}
}
//***************************************************************************************************************
int ChimeraCheckRDP::getChimeras(Sequence* query) {
	try {
		
		IS.clear();
				
		querySeq = query;
			
		closest = templateDB->findClosestSequence(query);  
	
		IS = findIS(); 
					
		//determine chimera report cutoff - window score above 95%
		//getCutoff();  - not very acurate predictor
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
vector<sim> ChimeraCheckRDP::findIS() {
	try {
		
		
		vector< map<int, int> > queryKmerInfo;	//vector of maps - each entry in the vector is a map of the kmers up to that spot in the unaligned seq
												//example:  seqKmerInfo[50] = map containing the kmers found in the first 50 + kmersize characters of ecoli.
												//i chose to store the kmers numbers in a map so you wouldn't have to check for dupilcate entries and could easily find the 
												//kmers 2 seqs had in common.  There may be a better way to do this thats why I am leaving so many comments...
		vector< map<int, int> > subjectKmerInfo;
		
		vector<sim>  isValues;
		string queryName = querySeq->getName();
		string seq = querySeq->getUnaligned();
		
		queryKmerInfo = kmer->getKmerCounts(seq);
		subjectKmerInfo = kmer->getKmerCounts(closest.getUnaligned());
		
		//find total kmers you have in common with closest[query] by looking at the last entry in the vector of maps for each
		int nTotal = calcKmers(queryKmerInfo[(queryKmerInfo.size()-1)], subjectKmerInfo[(subjectKmerInfo.size()-1)]);

		//you don't want the starting point to be virtually at the end so move it in 10%
		int start = seq.length() / 10;
        
		//for each window
		for (int f = start; f < (seq.length() - start); f+=increment) {
		
            if ((f - kmerSize) < 0)  { m->mothurOut("[ERROR]: Sequence " + querySeq->getName() + " is too short for your kmerSize, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
			
            if (m->control_pressed) { return isValues; }
            
			sim temp;
			
			string fragLeft = seq.substr(0, f);  //left side of breakpoint
			string fragRight = seq.substr(f);  //right side of breakpoint
			
			//make a sequence of the left side and right side
			Sequence* left = new Sequence(queryName, fragLeft);
			Sequence* right = new Sequence(queryName, fragRight);
			
			//find seqs closest to each fragment
			Sequence closestLeft = templateDB->findClosestSequence(left); 
	
			Sequence closestRight = templateDB->findClosestSequence(right); 
		
			//get kmerinfo for the closest left
			vector< map<int, int> > closeLeftKmerInfo = kmer->getKmerCounts(closestLeft.getUnaligned());
			
			//get kmerinfo for the closest right
			vector< map<int, int> > closeRightKmerInfo = kmer->getKmerCounts(closestRight.getUnaligned());
			
			//right side is tricky - since the counts grow on eachother to find the correct counts of only the right side you must subtract the counts of the left side
			//iterate through left sides map to subtract the number of times you saw things before you got the the right side
			map<int, int> rightside = queryKmerInfo[queryKmerInfo.size()-1];
			for (map<int, int>::iterator itleft = queryKmerInfo[f-kmerSize].begin(); itleft != queryKmerInfo[f-kmerSize].end(); itleft++) {
				int howManyTotal = queryKmerInfo[queryKmerInfo.size()-1][itleft->first];   //times that kmer was seen in total

				//itleft->second is times it was seen in left side, so howmanytotal - leftside should give you right side
				int howmanyright = howManyTotal - itleft->second;
				
				//if any were seen just on the left erase
				if (howmanyright == 0) {
					rightside.erase(itleft->first);
				}
			}
			
			map<int, int> closerightside = closeRightKmerInfo[closeRightKmerInfo.size()-1];
			for (map<int, int>::iterator itright = closeRightKmerInfo[f-kmerSize].begin(); itright != closeRightKmerInfo[f-kmerSize].end(); itright++) {
				int howManyTotal = closeRightKmerInfo[(closeRightKmerInfo.size()-1)][itright->first];   //times that kmer was seen in total

				//itleft->second is times it was seen in left side, so howmanytotal - leftside should give you right side
				int howmanyright = howManyTotal - itright->second;
				
				//if any were seen just on the left erase
				if (howmanyright == 0) {
					closerightside.erase(itright->first);
				}
			}

			
			int nLeft = calcKmers(closeLeftKmerInfo[f-kmerSize], queryKmerInfo[f-kmerSize]);

			int nRight = calcKmers(closerightside, rightside);

			int is = nLeft + nRight - nTotal;

			//save IS, leftparent, rightparent, breakpoint
			temp.leftParent = closestLeft.getName();
			temp.rightParent = closestRight.getName();
			temp.score = is;
			temp.midpoint = f;
			
			isValues.push_back(temp);
			
			delete left;
			delete right;
		}
		
		return isValues;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "findIS");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraCheckRDP::readName(string namefile) {
	try{
	
		string name;

		ifstream in;
		m->openInputFile(namefile, in);
				
		while (!in.eof()) {
			in >> name; m->gobble(in);
			names[name] = name;
		}
		in.close();
	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "readName");
		exit(1);
	}
}

//***************************************************************************************************************
//find the smaller map and iterate through it and count kmers in common
int ChimeraCheckRDP::calcKmers(map<int, int> query, map<int, int> subject) {
	try{
		
		int common = 0;
		
		map<int, int>::iterator smallone;
		map<int, int>::iterator largeone;

		if (query.size() < subject.size()) {
		
			for (smallone = query.begin(); smallone != query.end(); smallone++) {
				largeone = subject.find(smallone->first);
				
				//if you found it they have that kmer in common
				if (largeone != subject.end()) {	common++;	}
			}
			
		}else { 
		 
			for (smallone = subject.begin(); smallone != subject.end(); smallone++) {
				largeone = query.find(smallone->first);
				
				//if you found it they have that kmer in common
				if (largeone != query.end()) {		common++;	 }
			}
		}
		
		return common;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "calcKmers");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraCheckRDP::makeSVGpic(vector<sim> info) {
	try{
		
		string file = outputDir + querySeq->getName() + ".chimeracheck.svg";
		ofstream outsvg;
		m->openOutputFile(file, outsvg);
		
		int width = (info.size()*5) + 150;
		
		outsvg << "<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\" viewBox=\"0 0 700 " + toString(width) + "\">\n";
		outsvg << "<g>\n";
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString((width / 2) - 150) + "\" y=\"25\">Plotted IS values for " + querySeq->getName() + "</text>\n";
		
		outsvg <<  "<line x1=\"75\" y1=\"600\" x2=\"" + toString((info.size()*5) + 75) + "\" y2=\"600\" stroke=\"black\" stroke-width=\"2\"/>\n";  
		outsvg <<  "<line x1=\"75\" y1=\"600\" x2=\"75\" y2=\"125\" stroke=\"black\" stroke-width=\"2\"/>\n";
		
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"80\" y=\"620\">" + toString(info[0].midpoint) + "</text>\n";
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString((info.size()*5) + 75) + "\" y=\"620\">" + toString(info[info.size()-1].midpoint) + "</text>\n";
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"" + toString((width / 2) - 150) + "\" y=\"650\">Base Positions</text>\n";
		
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"50\" y=\"580\">0</text>\n";
		
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"50\" y=\"350\">IS</text>\n";
		
		
		//find max is score
		float biggest = 0.0;
		for (int i = 0; i < info.size(); i++) {
			if (info[i].score > biggest)  {
				biggest = info[i].score;
			}
		}
		
		outsvg << "<text fill=\"black\" class=\"seri\" x=\"50\" y=\"135\">" + toString(biggest) + "</text>\n";
		
		int scaler2 = 500 / biggest;
		
		
		outsvg << "<polyline fill=\"none\" stroke=\"red\" stroke-width=\"2\" points=\"";
		//160,200 180,230 200,210 234,220\"/> "; 
		for (int i = 0; i < info.size(); i++) {
			if(info[i].score < 0) { info[i].score = 0; }
			outsvg << ((i*5) + 75) << "," << (600 - (info[i].score * scaler2)) << " ";
		}
		
		outsvg << "\"/> ";
		outsvg << "</g>\n</svg>\n";
		
		outsvg.close();

	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckRDP", "makeSVGpic");
		exit(1);
	}
}
//***************************************************************************************************************/


