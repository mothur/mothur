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
#include "eachgapdist.h"

//********************************************************************************************************************
//sorts lowest to highest
inline bool compareQuanMembers(quanMember left, quanMember right){
	return (left.score < right.score);	
} 
//***************************************************************************************************************

Pintail::Pintail(string filename, string temp, bool f, int p, string mask, string cons, string q, int win, int inc, string o) : MothurChimera() { 
	try {
	
		fastafile = filename; 
		templateFileName = temp; templateSeqs = readSeqs(temp);
		filter = f;
		processors = p;
		setMask(mask);
		consfile = cons;
		quanfile = q;
		window = win;
		increment = inc; 
		outputDir = o; 
		
		distcalculator = new eachGapDist();
		decalc = new DeCalculator();
		
		doPrep();
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "Pintail");
		exit(1);
	}

}
//***************************************************************************************************************

Pintail::~Pintail() {
	try {
		
		delete distcalculator;
		delete decalc; 
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "~Pintail");
		exit(1);
	}
}
//***************************************************************************************************************
int Pintail::doPrep() {
	try {
		
		mergedFilterString = "";
		windowSizesTemplate.resize(templateSeqs.size(), window);
		quantiles.resize(100);  //one for every percent mismatch
		quantilesMembers.resize(100);  //one for every percent mismatch
		
		//if the user does not enter a mask then you want to keep all the spots in the alignment
		if (seqMask.length() == 0)	{	decalc->setAlignmentLength(templateSeqs[0]->getAligned().length());	}
		else						{	decalc->setAlignmentLength(seqMask.length());						}
		
		decalc->setMask(seqMask);

		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
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
			templateLines.push_back(new linePair(0, templateSeqs.size()));
		#endif
		
		m->mothurOut("Getting conservation... "); cout.flush();
		if (consfile == "") { 
			m->mothurOut("Calculating probability of conservation for your template sequences.  This can take a while...  I will output the frequency of the highest base in each position to a .freq file so that you can input them using the conservation parameter next time you run this command.  Providing the .freq file will improve speed.    "); cout.flush();
			probabilityProfile = decalc->calcFreq(templateSeqs, templateFileName); 
			if (m->getControl_pressed()) {  return 0;  }
			m->mothurOut("Done."); m->mothurOutEndLine();
		}else				{   probabilityProfile = readFreq();	m->mothurOut("Done.");		  }
		m->mothurOutEndLine();
		
		//make P into Q
		for (int i = 0; i < probabilityProfile.size(); i++)  { probabilityProfile[i] = 1 - probabilityProfile[i];  }  //
		
		bool reRead = false;
		//create filter if needed for later
		if (filter) {
						
			//read in all query seqs
			vector<Sequence*> tempQuerySeqs = readSeqs(fastafile);
				
			vector<Sequence*> temp;
			//merge query seqs and template seqs
			temp = templateSeqs;
			for (int i = 0; i < tempQuerySeqs.size(); i++) {  temp.push_back(tempQuerySeqs[i]);  }
	
			if (seqMask != "") {
			    reRead = true;
				//mask templates
				for (int i = 0; i < temp.size(); i++) {
					if (m->getControl_pressed()) {  
						for (int i = 0; i < tempQuerySeqs.size(); i++) { delete tempQuerySeqs[i];  }
						return 0; 
					}
					decalc->runMask(temp[i]);
				}
			}

			mergedFilterString = createFilter(temp, 0.5);
			
			if (m->getControl_pressed()) {  
				for (int i = 0; i < tempQuerySeqs.size(); i++) { delete tempQuerySeqs[i];  }
				return 0; 
			}
			
			//reread template seqs
			for (int i = 0; i < tempQuerySeqs.size(); i++) { delete tempQuerySeqs[i];  }
		}
		
		
		//quantiles are used to determine whether the de values found indicate a chimera
		//if you have to calculate them, its time intensive because you are finding the de and deviation values for each 
		//combination of sequences in the template
		if (quanfile != "") {  
			quantiles = readQuantiles(); 
		}else {
			if ((!filter) && (seqMask != "")) { //if you didn't filter but you want to mask. if you filtered then you did mask first above.
				reRead = true;
				//mask templates
				for (int i = 0; i < templateSeqs.size(); i++) {
					if (m->getControl_pressed()) {  return 0;  }
					decalc->runMask(templateSeqs[i]);
				}
			}
			
			if (filter) { 
				reRead = true;
				for (int i = 0; i < templateSeqs.size(); i++) {
					if (m->getControl_pressed()) {  return 0;  }
					runFilter(templateSeqs[i]);
				}
			}
			
			m->mothurOut("Calculating quantiles for your template.  This can take a while...  I will output the quantiles to a .quan file that you can input them using the quantiles parameter next time you run this command.  Providing the .quan file will dramatically improve speed.    "); cout.flush();
			if (processors == 1) { 
				quantilesMembers = decalc->getQuantiles(templateSeqs, windowSizesTemplate, window, probabilityProfile, increment, 0, templateSeqs.size());
			}else {		createProcessesQuan();		}
		
			if (m->getControl_pressed()) {  return 0;  }
			
			string noOutliers, outliers;
			
			if ((!filter) && (seqMask == "")) {
				noOutliers = m->getRootName(m->getSimpleName(templateFileName)) + "pintail.quan";
			}else if ((!filter) && (seqMask != "")) { 
				noOutliers =m->getRootName(m->getSimpleName(templateFileName)) + "pintail.masked.quan";
			}else if ((filter) && (seqMask != "")) { 
				noOutliers = m->getRootName(m->getSimpleName(templateFileName)) + "pintail.filtered." + m->getSimpleName(m->getRootName(fastafile)) + "masked.quan";
			}else if ((filter) && (seqMask == "")) { 
				noOutliers = m->getRootName(m->getSimpleName(templateFileName)) + "pintail.filtered." + m->getSimpleName(m->getRootName(fastafile)) + "quan";
			}

			decalc->removeObviousOutliers(quantilesMembers, templateSeqs.size());
			
			if (m->getControl_pressed()) {  return 0;  }
		
			string outputString = "#" + m->getVersion() + "\n";
			
			//adjust quantiles
			for (int i = 0; i < quantilesMembers.size(); i++) {
				vector<float> temp;
				
				if (quantilesMembers[i].size() == 0) {
					//in case this is not a distance found in your template files
					for (int g = 0; g < 6; g++) {
						temp.push_back(0.0);
					}
				}else{
					
					sort(quantilesMembers[i].begin(), quantilesMembers[i].end());
					
					//save 10%
					temp.push_back(quantilesMembers[i][int(quantilesMembers[i].size() * 0.10)]);
					//save 25%
					temp.push_back(quantilesMembers[i][int(quantilesMembers[i].size() * 0.25)]);
					//save 50%
					temp.push_back(quantilesMembers[i][int(quantilesMembers[i].size() * 0.5)]);
					//save 75%
					temp.push_back(quantilesMembers[i][int(quantilesMembers[i].size() * 0.75)]);
					//save 95%
					temp.push_back(quantilesMembers[i][int(quantilesMembers[i].size() * 0.95)]);
					//save 99%
					temp.push_back(quantilesMembers[i][int(quantilesMembers[i].size() * 0.99)]);
					
				}
				
				//output quan value
				outputString += toString(i+1);
				for (int u = 0; u < temp.size(); u++) {   outputString += "\t" + toString(temp[u]); }
				outputString += "\n";
				
				quantiles[i] = temp;
				
			}
			
			printQuanFile(noOutliers, outputString);
			
			//free memory
			quantilesMembers.clear();
			
			m->mothurOut("Done."); m->mothurOutEndLine();
		}
		
		if (reRead) {
			for (int i = 0; i < templateSeqs.size(); i++) { delete templateSeqs[i];  }
			templateSeqs.clear();
			templateSeqs = readSeqs(templateFileName);
		}

		
		//free memory
		for (int i = 0; i < templateLines.size(); i++) { delete templateLines[i];  }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "doPrep");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence Pintail::print(ostream& out, ostream& outAcc) {
	try {
		
		int index = ceil(deviation);
		
		//is your DE value higher than the 95%
		string chimera;
		if (index != 0) {  //if index is 0 then its an exact match to a template seq
			if (quantiles[index][4] == 0.0) {
				chimera = "Your template does not include sequences that provide quantile values at distance " + toString(index);
			}else {
				if (DE > quantiles[index][4])		{	chimera = "Yes";	}
				else								{	chimera = "No";		}
			}
		}else{ chimera = "No";		}
		
		out << querySeq->getName() << '\t' << "div: " << deviation << "\tstDev: " << DE << "\tchimera flag: " << chimera << endl;
		if (chimera == "Yes") {
			m->mothurOut(querySeq->getName() + "\tdiv: " + toString(deviation) + "\tstDev: " + toString(DE) + "\tchimera flag: " + chimera); m->mothurOutEndLine();
			outAcc << querySeq->getName() << endl;
		}
		out << "Observed";
		
		for (int j = 0; j < obsDistance.size(); j++) {  out  << '\t' << obsDistance[j];  }
		out << endl;
		
		out << "Expected";
		
		for (int m = 0; m < expectedDistance.size(); m++) {  out << '\t' << expectedDistance[m] ;  }
		out << endl;
		
		return *querySeq;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "print");
		exit(1);
	}
}
//***************************************************************************************************************
int Pintail::getChimeras(Sequence* query) {
	try {
		querySeq = query;
		trimmed.clear();
		windowSizes = window;
							
		//find pairs has to be done before a mask
		bestfit = findPairs(query);
		
		if (m->getControl_pressed()) {  return 0; } 
		
		//if they mask  
		if (seqMask != "") {
			decalc->runMask(query);
			decalc->runMask(bestfit);
		}

		if (filter) { //must be done after a mask
			runFilter(query);
			runFilter(bestfit);
		}
		
				
		//trim seq
		decalc->trimSeqs(query, bestfit, trimmed);  
		
		//find windows
		it = trimmed.begin();
		windowsForeachQuery = decalc->findWindows(query, it->first, it->second, windowSizes, increment);

		//find observed distance
		obsDistance = decalc->calcObserved(query, bestfit, windowsForeachQuery, windowSizes);
		
		if (m->getControl_pressed()) {  return 0; } 
				
		Qav = decalc->findQav(windowsForeachQuery, windowSizes, probabilityProfile);
		
		if (m->getControl_pressed()) {  return 0; } 

		//find alpha			
		seqCoef	= decalc->getCoef(obsDistance, Qav);
		
		//calculating expected distance
		expectedDistance = decalc->calcExpected(Qav, seqCoef);
		
		if (m->getControl_pressed()) {  return 0; } 
		
		//finding de
		DE = decalc->calcDE(obsDistance, expectedDistance);
		
		if (m->getControl_pressed()) {  return 0; } 
		
		//find distance between query and closest match
		it = trimmed.begin();
		deviation = decalc->calcDist(query, bestfit, it->first, it->second); 
		
		delete bestfit;
									
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "getChimeras");
		exit(1);
	}
}

//***************************************************************************************************************

vector<float> Pintail::readFreq() {
	try {
		//read in probabilities and store in vector
		int pos; float num; 

		vector<float> prob;
		set<int> h = decalc->getPos();  //positions of bases in masking sequence
		
		ifstream in;
		m->openInputFile(consfile, in);
		
		//read version
		string line = m->getline(in); m->gobble(in);
				
		while(!in.eof()){
			
			in >> pos >> num;
			
			if (h.count(pos) > 0) {
				float Pi;
				Pi =  (num - 0.25) / 0.75; 
			
				//cannot have probability less than 0.
				if (Pi < 0) { Pi = 0.0; }

				//do you want this spot
				prob.push_back(Pi);  
			}
			
			m->gobble(in);
		}
		in.close();
		return prob;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "readFreq");
		exit(1);
	}
}

//***************************************************************************************************************
//calculate the distances from each query sequence to all sequences in the template to find the closest sequence
Sequence* Pintail::findPairs(Sequence* q) {
	try {
		
		Sequence* seqsMatches;  
		
		seqsMatches = decalc->findClosest(q, templateSeqs);
		return seqsMatches;
	
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "findPairs");
		exit(1);
	}
}
//**************************************************************************************************
void Pintail::createProcessesQuan() {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 1;
		vector<int> processIDS;
        bool recalc = false;
				
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  
				process++;
			}else if (pid == 0){
				
				quantilesMembers = decalc->getQuantiles(templateSeqs, windowSizesTemplate, window, probabilityProfile, increment, templateLines[process]->start, templateLines[process]->end);
				
				//write out data to file so parent can read it
				ofstream out;
				string s = m->mothurGetpid(process) + ".temp";
				m->openOutputFile(s, out);
								
				//output observed distances
				for (int i = 0; i < quantilesMembers.size(); i++) {
					out << quantilesMembers[i].size();
					for (int j = 0; j < quantilesMembers[i].size(); j++) {
						out << '\t' << quantilesMembers[i][j];
					}
					out << endl;
				}
				
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->setControl_pressed(false);
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove((toString(processIDS[i]) + ".temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->setControl_pressed(false);
					 for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            for (int i = 0; i < templateLines.size(); i++) {  delete templateLines[i];  }  templateLines.clear();
            for (int i = 0; i < processors; i++) {
                templateLines.push_back(new linePair());
                templateLines[i]->start = int (sqrt(float(i)/float(processors)) * templateSeqs.size());
                templateLines[i]->end = int (sqrt(float(i+1)/float(processors)) * templateSeqs.size());
            }
            
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);
                    process++;
                }else if (pid == 0){
                    
                    quantilesMembers = decalc->getQuantiles(templateSeqs, windowSizesTemplate, window, probabilityProfile, increment, templateLines[process]->start, templateLines[process]->end);
                    
                    //write out data to file so parent can read it
                    ofstream out;
                    string s = m->mothurGetpid(process) + ".temp";
                    m->openOutputFile(s, out);
                    
                    //output observed distances
                    for (int i = 0; i < quantilesMembers.size(); i++) {
                        out << quantilesMembers[i].size();
                        for (int j = 0; j < quantilesMembers[i].size(); j++) {
                            out << '\t' << quantilesMembers[i][j];
                        }
                        out << endl;
                    }
                    
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

		
		//parent does its part
		quantilesMembers = decalc->getQuantiles(templateSeqs, windowSizesTemplate, window, probabilityProfile, increment, templateLines[0]->start, templateLines[0]->end);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<(processors-1);i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}

		//get data created by processes
		for (int i=0;i<(processors-1);i++) { 
			ifstream in;
			string s = toString(processIDS[i]) + ".temp";
			m->openInputFile(s, in);
			
			vector< vector<float> > quan; 
			quan.resize(100);
			
			//get quantiles
			for (int h = 0; h < quan.size(); h++) {
				int num;
				in >> num; 
				
				m->gobble(in);

				vector<float> q;  float w; 
				for (int j = 0; j < num; j++) {
					in >> w;
					q.push_back(w);
				}

				quan[h] = q;
				m->gobble(in);
			}
			
	
			//save quan in quantiles
			for (int j = 0; j < quan.size(); j++) {
				//put all values of q[i] into quan[i]
				for (int l = 0; l < quan[j].size(); l++) {  quantilesMembers[j].push_back(quan[j][l]);   }
				//quantilesMembers[j].insert(quantilesMembers[j].begin(), quan[j].begin(), quan[j].end());
			}
					
			in.close();
			m->mothurRemove(s);
		}

#else
		quantilesMembers = decalc->getQuantiles(templateSeqs, windowSizesTemplate, window, probabilityProfile, increment, 0, templateSeqs.size());
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "createProcessesQuan");
		exit(1);
	}
}
//***************************************************************************************************************
vector< vector<float> > Pintail::readQuantiles() {
	try {
		int num; 
		float ten, twentyfive, fifty, seventyfive, ninetyfive, ninetynine; 
		
		vector< vector<float> > quan;
		vector <float> temp; temp.resize(6, 0);
		
		//to fill 0
		quan.push_back(temp); 

		ifstream in;
		m->openInputFile(quanfile, in);
		
		//read version
		string line = m->getline(in); m->gobble(in);
			
		while(!in.eof()){
			
			in >> num >> ten >> twentyfive >> fifty >> seventyfive >> ninetyfive >> ninetynine; 
			
			temp.clear();
			
			temp.push_back(ten); 
			temp.push_back(twentyfive);
			temp.push_back(fifty);
			temp.push_back(seventyfive);
			temp.push_back(ninetyfive);
			temp.push_back(ninetynine);
			
			quan.push_back(temp);  
	
			m->gobble(in);
		}
		in.close();

		return quan;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "readQuantiles");
		exit(1);
	}
}
//***************************************************************************************************************/

void Pintail::printQuanFile(string file, string outputString) {
	try {
	
			ofstream outQuan;
			m->openOutputFile(file, outQuan);
			
			outQuan << outputString;
			
			outQuan.close();
	}
	catch(exception& e) {
		m->errorOut(e, "Pintail", "printQuanFile");
		exit(1);
	}
}

//***************************************************************************************************************/



