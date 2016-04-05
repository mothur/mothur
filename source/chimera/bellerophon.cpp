/*
 *  bellerophon.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "bellerophon.h"
#include "eachgapdist.h"
#include "ignoregaps.h"
#include "onegapdist.h"


/***************************************************************************************************************/

Bellerophon::Bellerophon(string name, bool filterSeqs,  bool c, int win, int inc, int p, string o) : MothurChimera() {
	try {
		fastafile = name;
		correction = c;
		outputDir = o;
		window = win;
		increment = inc;
		processors = p;
		
		//read in sequences
		seqs = readSeqs(fastafile);
		numSeqs = seqs.size();
		if (numSeqs == 0) { m->mothurOut("Error in reading you sequences."); m->mothurOutEndLine(); exit(1); }
	
		//do soft filter
		if (filterSeqs)  {
			createFilter(seqs, 0.5);
			for (int i = 0; i < seqs.size(); i++) {  runFilter(seqs[i]);  }
		}
		
		distCalculator = new eachGapDist();
		
		//set default window to 25% of sequence length
		string seq0 = seqs[0]->getAligned();
		if (window == 0) { window = seq0.length() / 4;  }
		else if (window > (seq0.length() / 2)) {  
			m->mothurOut("Your sequence length is = " + toString(seq0.length()) + ". You have selected a window size greater than the length of half your aligned sequence. I will run it with a window size of " + toString((seq0.length() / 2))); m->mothurOutEndLine();
			window = (seq0.length() / 2);
		}
		
		if (increment > (seqs[0]->getAlignLength() - (2*window))) { 
			if (increment != 10) {
			
				m->mothurOut("You have selected a increment that is too large. I will use the default."); m->mothurOutEndLine();
				increment = 10;
				if (increment > (seqs[0]->getAlignLength() - (2*window))) {  increment = 0;  }
				
			}else{ increment = 0; }
		}
		
		if (increment == 0) { iters = 1; }
		else { iters = ((seqs[0]->getAlignLength() - (2*window)) / increment); }
		
		//initialize pref
		pref.resize(iters);
		for (int i = 0; i < iters; i++) { 
			Preference temp;
			for (int j = 0; j < numSeqs; j++) {  
				pref[i].push_back(temp); 
			}
		} 

	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "Bellerophon");
		exit(1);
	}
}

//***************************************************************************************************************
int Bellerophon::print(ostream& out, ostream& outAcc, string s) {
	try {
		int above1 = 0;
		
		//sorted "best" preference scores for all seqs
		vector<Preference> best = getBestPref();
		
		if (m->control_pressed) { return numSeqs; }
		
		out << "Name\tScore\tLeft\tRight\t" << endl;
		//output prefenence structure to .chimeras file
		for (int i = 0; i < best.size(); i++) {
			
			if (m->control_pressed) {  return numSeqs; }
			
			out << best[i].name << '\t' << setprecision(3) << best[i].score << '\t' << best[i].leftParent << '\t' << best[i].rightParent << endl;
			
			//calc # of seqs with preference above 1.0
			if (best[i].score > 1.0) { 
				above1++; 
				outAcc << best[i].name << endl;
				m->mothurOut(best[i].name + " is a suspected chimera at breakpoint " + toString(best[i].midpoint)); m->mothurOutEndLine();
				m->mothurOut("It's score is " + toString(best[i].score) + " with suspected left parent " + best[i].leftParent + " and right parent " + best[i].rightParent); m->mothurOutEndLine();
			}
		}
		
		//output results to screen
		m->mothurOutEndLine();
		m->mothurOut("Sequence with preference score above 1.0: " + toString(above1)); m->mothurOutEndLine();
		int spot;
		spot = best.size()-1;
		m->mothurOut("Minimum:\t" + toString(best[spot].score)); m->mothurOutEndLine();
		spot = best.size() * 0.975;
		m->mothurOut("2.5%-tile:\t" + toString(best[spot].score)); m->mothurOutEndLine();
		spot = best.size() * 0.75;
		m->mothurOut("25%-tile:\t" + toString(best[spot].score)); m->mothurOutEndLine();
		spot = best.size() * 0.50;
		m->mothurOut("Median: \t" + toString(best[spot].score)); m->mothurOutEndLine();
		spot = best.size() * 0.25;
		m->mothurOut("75%-tile:\t" + toString(best[spot].score)); m->mothurOutEndLine();
		spot = best.size() * 0.025;
		m->mothurOut("97.5%-tile:\t" + toString(best[spot].score)); m->mothurOutEndLine();
		spot = 0;
		m->mothurOut("Maximum:\t" + toString(best[spot].score)); m->mothurOutEndLine();
		
		return numSeqs;

	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "print");
		exit(1);
	}
}
//********************************************************************************************************************
//sorts highest score to lowest
inline bool comparePref(Preference left, Preference right){
	return (left.score > right.score);	
}
//***************************************************************************************************************
int Bellerophon::getChimeras() {
	try {
		
		//create breaking points
		vector<int> midpoints;   midpoints.resize(iters, window);
		for (int i = 1; i < iters; i++) {  midpoints[i] = midpoints[i-1] + increment;  }
	
		
		//divide breakpoints between processors
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			if(processors == 1){ 
				lines.push_back(linePair(0, iters));	
				
				//fill pref with scores
				driverChimeras(midpoints, lines[0]);
	
			}else{
			
				int numSeqsPerProcessor = iters / processors;
				
				for (int i = 0; i < processors; i++) {
					unsigned long long startPos = i * numSeqsPerProcessor;
					if(i == processors - 1){
						numSeqsPerProcessor = iters - i * numSeqsPerProcessor;
					}
					lines.push_back(linePair(startPos, numSeqsPerProcessor));
				}
				
				createProcesses(midpoints);
			}
		#else
			lines.push_back(linePair(0, iters));	
			
			///fill pref with scores
			driverChimeras(midpoints, lines[0]);
		#endif
	
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "getChimeras");
		exit(1);
	}
}
/**************************************************************************************************/

int Bellerophon::createProcesses(vector<int> mid) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 0;
		int exitCommand = 1;
		vector<int> processIDS;
        bool recalc = false;
				
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				exitCommand = driverChimeras(mid, lines[process]);
				string tempOut = outputDir + toString(m->mothurGetpid(process)) + ".temp";
				writePrefs(tempOut, lines[process]);
				exit(0);
            }else {
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                m->control_pressed = false;
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                recalc = true;
                break;
            }
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); }  m->control_pressed = false; processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            lines.clear();
            int numSeqsPerProcessor = iters / processors;
            
            for (int i = 0; i < processors; i++) {
                unsigned long long startPos = i * numSeqsPerProcessor;
                if(i == processors - 1){
                    numSeqsPerProcessor = iters - i * numSeqsPerProcessor;
                }
                lines.push_back(linePair(startPos, numSeqsPerProcessor));
            }

            processIDS.clear();
            process = 0;
            
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    exitCommand = driverChimeras(mid, lines[process]);
                    string tempOut = outputDir + toString(m->mothurGetpid(process)) + ".temp";
                    writePrefs(tempOut, lines[process]);
                    exit(0);
                }else {
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) {
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//get info that other processes created
		for (int i = 0; i < processIDS.size(); i++) {
			string tempIn = outputDir + toString(processIDS[i]) + ".temp";
			readPrefs(tempIn);
		}
		
		return exitCommand;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "createProcesses");
		exit(1);
	}
}
//***************************************************************************************************************
int Bellerophon::driverChimeras(vector<int> midpoints, linePair line) {
	try {
		
		for (int h = line.start; h < (line.start + line.num); h++) {
			count = h;
			int midpoint = midpoints[h];
		
			//initialize pref[count]		
			for (int i = 0; i < numSeqs; i++ ) { 
				pref[count][i].name = seqs[i]->getName();
				pref[count][i].midpoint = midpoint;  
			}
			
			if (m->control_pressed) { return 0; }
			
			//create 2 vectors of sequences, 1 for left side and one for right side
			vector<Sequence> left;  vector<Sequence> right;
			
			for (int i = 0; i < seqs.size(); i++) {
				
				if (m->control_pressed) { return 0; }
				
				//cout << "midpoint = " << midpoint << "\twindow = " << window << endl;
				//cout << "whole = " << seqs[i]->getAligned().length() << endl;
				//save left side
				string seqLeft = seqs[i]->getAligned().substr(midpoint-window, window);
				Sequence tempLeft;
				tempLeft.setName(seqs[i]->getName());
				tempLeft.setAligned(seqLeft);
				left.push_back(tempLeft);
				//cout << "left = " << tempLeft.getAligned().length() << endl;			
				//save right side
				string seqRight = seqs[i]->getAligned().substr(midpoint, window);
				Sequence tempRight;
				tempRight.setName(seqs[i]->getName());
				tempRight.setAligned(seqRight);
				right.push_back(tempRight);
				//cout << "right = " << seqRight.length() << endl;	
			}
			
			//this should be parallelized
			//perference = sum of (| distance of my left to sequence j's left - distance of my right to sequence j's right | )
			//create a matrix containing the distance from left to left and right to right
			//calculate distances
			SparseMatrix* SparseLeft = new SparseMatrix();
			SparseMatrix* SparseRight = new SparseMatrix();
			
			createSparseMatrix(0, left.size(), SparseLeft, left);
			
			if (m->control_pressed) { delete SparseLeft; delete SparseRight; return 0; }
			
			createSparseMatrix(0, right.size(), SparseRight, right);
			
			if (m->control_pressed) { delete SparseLeft; delete SparseRight; return 0; }
			
			left.clear(); right.clear();
			vector<SeqMap> distMapRight;
			vector<SeqMap> distMapLeft;
			
			// Create a data structure to quickly access the distance information.
			//this is from thallingers reimplementation on get.oturep
			// It consists of a vector of distance maps, where each map contains
			// all distances of a certain sequence. Vector and maps are accessed
			// via the index of a sequence in the distance matrix
			distMapRight = vector<SeqMap>(numSeqs); 
			distMapLeft = vector<SeqMap>(numSeqs); 
			//cout << "left" << endl << endl;
			for (MatData currentCell = SparseLeft->begin(); currentCell != SparseLeft->end(); currentCell++) {
				distMapLeft[currentCell->row][currentCell->column] = currentCell->dist;
				if (m->control_pressed) { delete SparseLeft; delete SparseRight; return 0; }
				//cout << " i = " << currentCell->row << " j = " << currentCell->column << " dist = " << currentCell->dist << endl;
			}
			//cout << "right" << endl << endl;
			for (MatData currentCell = SparseRight->begin(); currentCell != SparseRight->end(); currentCell++) {
				distMapRight[currentCell->row][currentCell->column] = currentCell->dist;
				if (m->control_pressed) { delete SparseLeft; delete SparseRight; return 0; }
				//cout << " i = " << currentCell->row << " j = " << currentCell->column << " dist = " << currentCell->dist << endl;
			}
			
			delete SparseLeft;
			delete SparseRight;
			
			//fill preference structure
			generatePreferences(distMapLeft, distMapRight, midpoint);
			
			if (m->control_pressed) { return 0; }
			
			//report progress
			if((h+1) % 10 == 0){ m->mothurOutJustToScreen("Processing sliding window: " + toString(h+1) + "\n") ;		}
			
		}
		
		//report progress
		if((line.start + line.num) % 10 != 0){  m->mothurOutJustToScreen("Processing sliding window: " + toString(line.start + line.num) + "\n") ;		}

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "driverChimeras");
		exit(1);
	}
}

/***************************************************************************************************************/
int Bellerophon::createSparseMatrix(int startSeq, int endSeq, SparseMatrix* sparse, vector<Sequence> s){
	try {

		for(int i=startSeq; i<endSeq; i++){
			
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) { return 0; }
			
				distCalculator->calcDist(s[i], s[j]);
				float dist = distCalculator->getDist();
			
				PCell temp(i, j, dist);
				sparse->addCell(temp);
				
			}
		}
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "createSparseMatrix");
		exit(1);
	}
}
/***************************************************************************************************************/
int Bellerophon::generatePreferences(vector<SeqMap> left, vector<SeqMap> right, int mid){
	try {
		
		SeqMap::iterator itR;
		SeqMap::iterator itL;
		
		for (int i = 0; i < left.size(); i++) {
			
			SeqMap currentLeft = left[i];    //example i = 3;   currentLeft is a map of 0 to the distance of sequence 3 to sequence 0,
												//										1 to the distance of sequence 3 to sequence 1,
												//										2 to the distance of sequence 3 to sequence 2.
			SeqMap currentRight = right[i];		// same as left but with distances on the right side.
			
			for (int j = 0; j < i; j++) {
			
				if (m->control_pressed) {  return 0; }
				
				itL = currentLeft.find(j);
				itR = currentRight.find(j);
				
				//if you can find this entry update the preferences
				if ((itL != currentLeft.end()) && (itR != currentRight.end())) {
				
					if (!correction) {
						pref[count][i].score += abs((itL->second - itR->second));
						pref[count][j].score += abs((itL->second - itR->second));
					}else {
						pref[count][i].score += abs((sqrt(itL->second) - sqrt(itR->second)));
						pref[count][j].score += abs((sqrt(itL->second) - sqrt(itR->second)));
                    }
					//are you the closest left sequence
					if (itL->second < pref[count][i].closestLeft) {  

						pref[count][i].closestLeft = itL->second;
						pref[count][i].leftParent = seqs[j]->getName();
					}
					if (itL->second < pref[count][j].closestLeft) {
						pref[count][j].closestLeft = itL->second;
						pref[count][j].leftParent = seqs[i]->getName();
					}
					
					//are you the closest right sequence
					if (itR->second < pref[count][i].closestRight) {   
						pref[count][i].closestRight = itR->second;
						pref[count][i].rightParent = seqs[j]->getName();
					}
					if (itR->second < pref[count][j].closestRight) {   
						pref[count][j].closestRight = itR->second;
						pref[count][j].rightParent = seqs[i]->getName();
					}
					
				}
			}
		
		}
		
				
		return 1;

	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "generatePreferences");
		exit(1);
	}
}
/**************************************************************************************************/
vector<Preference> Bellerophon::getBestPref() {
	try {
		
		vector<Preference> best;
		
		//for each sequence
		for (int i = 0; i < numSeqs; i++) {
			
			//set best pref score to first one
			Preference temp = pref[0][i];
			
			if (m->control_pressed) { return best;  }
			
			//for each window
			for (int j = 1; j < pref.size(); j++) {
				
				//is this a better score
				if (pref[j][i].score > temp.score) {	temp = pref[j][i];		}
			}
			
			best.push_back(temp);
		}
		
		//rank preference score to eachother
		float dme = 0.0;
		float expectedPercent = 1 / (float) (best.size());
		
		for (int i = 0; i < best.size(); i++) {	 dme += best[i].score;  }
	
		for (int i = 0; i < best.size(); i++) {

			if (m->control_pressed) { return best; }
			
			//gives the actual percentage of the dme this seq adds
			best[i].score = best[i].score / dme;
			
			//how much higher or lower is this than expected
			best[i].score = best[i].score / expectedPercent;
		
		}
		
		//sort Preferences highest to lowest
		sort(best.begin(), best.end(), comparePref);

		return best;
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "getBestPref");
		exit(1);
	}
}
/**************************************************************************************************/
int Bellerophon::writePrefs(string file, linePair tempLine) {
	try {
	
		ofstream outTemp;
		m->openOutputFile(file, outTemp);
		
		//lets you know what part of the pref matrix you are writing
		outTemp << tempLine.start << '\t' << tempLine.num << endl;
		
		for (int i = tempLine.start; i < (tempLine.start + tempLine.num); i++) {
			
			for (int j = 0; j < numSeqs; j++) {
				
				if (m->control_pressed) { outTemp.close(); m->mothurRemove(file); return 0; }
				
				outTemp << pref[i][j].name << '\t' << pref[i][j].leftParent << '\t' << pref[i][j].rightParent << '\t';
				outTemp << pref[i][j].score << '\t' << pref[i][j].closestLeft << '\t' << pref[i][j].closestRight << '\t' << pref[i][j].midpoint <<  endl;
			}
		}
		
		outTemp.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "writePrefs");
		exit(1);
	}
}
/**************************************************************************************************/
int Bellerophon::readPrefs(string file) {
	try {
	
		ifstream inTemp;
		m->openInputFile(file, inTemp);
		
		int start, num;
		
		//lets you know what part of the pref matrix you are writing
		inTemp >> start >> num;  m->gobble(inTemp);
		
		for (int i = start; i < num; i++) {
			
			for (int j = 0; j < numSeqs; j++) {
				
				if (m->control_pressed) { inTemp.close(); m->mothurRemove(file); return 0; }
			
				inTemp >> pref[i][j].name >> pref[i][j].leftParent >> pref[i][j].rightParent;
				inTemp >> pref[i][j].score >> pref[i][j].closestLeft >> pref[i][j].closestRight >> pref[i][j].midpoint;
				m->gobble(inTemp);
			}
		}
		
		inTemp.close();
		
		m->mothurRemove(file);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "writePrefs");
		exit(1);
	}
}
/**************************************************************************************************/
vector<string> Bellerophon::getBestWindow(linePair line) {
	try {
	
		vector<string> best;
			
		//for each sequence
		for (int i = 0; i < numSeqs; i++) {
			
			//set best pref score to first one
			Preference temp = pref[line.start][i];
			
			if (m->control_pressed) { return best;  }
			
			//for each window
			for (int j = (line.start+1); j < (line.start+line.num); j++) {
				
				//is this a better score
				if (pref[j][i].score > temp.score) {	temp = pref[j][i];		}
			}
			
			string tempString = temp.name + '\t' + temp.leftParent + '\t' + temp.rightParent + '\t' + toString(temp.score);
			best.push_back(tempString);
		}

		return best;
	
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "getBestWindow");
		exit(1);
	}
}
/**************************************************************************************************/
int Bellerophon::fillPref(int process, vector<string>& best) {
	try {
		//figure out where you start so you can put the best scores there
		int numSeqsPerProcessor = iters / processors;
		int start = process * numSeqsPerProcessor;
		
		for (int i = 0; i < best.size(); i++) {
		
			if (m->control_pressed) { return 0;  }
			
			istringstream iss (best[i],istringstream::in);
			
			string tempScore;
			iss >> pref[start][i].name >> pref[start][i].leftParent >> pref[start][i].rightParent >> tempScore;
			convert(tempScore, pref[start][i].score); 
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Bellerophon", "fillPref");
		exit(1);
	}
}

/**************************************************************************************************/

