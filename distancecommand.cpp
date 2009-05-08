/*
 *  distancecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "distancecommand.h"
#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"

//**********************************************************************************************************************

DistanceCommand::DistanceCommand(){
	try {
		globaldata = GlobalData::getInstance();
		validCalculator = new ValidCalculators();
		ends = globaldata->getEnds();
		seqDB = globaldata->gSequenceDB;
		convert(globaldata->getProcessors(), processors);
		convert(globaldata->getCutOff(), cutoff);
		distFile = getRootName(globaldata->getFastaFile()) + "dist";
		
		int i;
		if (ends != "T") {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("distance", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "nogaps") { 
						distCalculator = new ignoreGaps();
					}else if (globaldata->Estimators[i] == "eachgap") { 
						distCalculator = new eachGapDist();	
					}else if (globaldata->Estimators[i] == "onegap") {
						distCalculator = new oneGapDist();					}
				}
			}
		}else {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("distance", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "nogaps") { 
						distCalculator = new ignoreGaps();	
					}else if (globaldata->Estimators[i] == "eachgap") { 
						distCalculator = new eachGapIgnoreTermGapDist();
					}else if (globaldata->Estimators[i] == "onegap") { 
						distCalculator = new oneGapIgnoreTermGapDist();	
					}
				}
			}
		}
				
				
		//reset calc for next command
		globaldata->setCalc("");

		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function DistanceCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function DistanceCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

int DistanceCommand::execute(){
	try {
		int numSeqs = seqDB->getNumSeqs();
		
		system(("rm "+distFile).c_str() );
		if(processors == 1){
			driver(distCalculator, seqDB, 0, numSeqs, distFile, cutoff);	
		}	
		else if(processors == 2){
		
			int pid = fork();
			if(pid > 0){
				driver(distCalculator, seqDB, 0, (numSeqs/sqrt(2)), distFile + "tempa", cutoff);	
//				system(("cat " + distFile + "tempa" + " >> " + distFile).c_str());
//				system(("rm " + distFile + "tempa").c_str());				
			}
			else{
				driver(distCalculator, seqDB, (numSeqs/sqrt(2)), numSeqs, distFile + "tempb", cutoff);	
//				system(("cat " + distFile + "tempb" + " >> " + distFile).c_str());
//				system(("rm " + distFile + "tempb").c_str());				
			}
			wait(NULL);

		}
		else if(processors == 3){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, sqrt(3) * numSeqs / 3, distFile + "tempa", cutoff);
					#ifdef HAVE_CAT
						system(("cat " + distFile + "tempa" + " >> " + distFile).c_str());
					#else
						#ifdef HAVE_COPY
//get system call from pat system(("copy " + distFile + "tempa").c_str());
						#else
							cout << "Sorry but I can't continue because this operating system doesn't appear to support the cat() or copy() system calls." << endl;
						#endif
					#endif
					
					#ifdef HAVE_RM
						system(("rm " + distFile + "tempa").c_str());	
					#else
						#ifdef HAVE_ERASE
							system(("erase " + distFile + "tempa").c_str());
						#else
							cout << "Sorry but I can't remove the required files because this operating system doesn't appear to support the rm() or erase() system calls." << endl;
						#endif
					#endif
				}
				else{
					driver(distCalculator, seqDB, sqrt(3) * numSeqs / 3, sqrt(6) * numSeqs / 3, distFile + "tempb", cutoff);	
					system(("cat " + distFile + "tempb" + " >> " + distFile).c_str());
					system(("rm " + distFile + "tempb").c_str());				
				}
				wait(NULL);
			}
			else{
				driver(distCalculator, seqDB, sqrt(6) * numSeqs / 3, numSeqs, distFile + "tempc", cutoff);	
				system(("cat " + distFile + "tempc" + " >> " + distFile).c_str());
				system(("rm " + distFile + "tempc").c_str());				
			}
			wait(NULL);
		}
		else if(processors == 4){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, numSeqs / 2, distFile + "tempa", cutoff);	
					system(("cat " + distFile + "tempa" + " >> " + distFile).c_str());
					system(("rm " + distFile + "tempa").c_str());				
				}
				else{
					driver(distCalculator, seqDB, numSeqs / 2, (numSeqs/sqrt(2)), distFile + "tempb", cutoff);	
					system(("cat " + distFile + "tempb" + " >> " + distFile).c_str());
					system(("rm " + distFile + "tempb").c_str());				
				}
				wait(NULL);
			}
			else{
				int pid3 = fork();
				if(pid3 > 0){
					driver(distCalculator, seqDB, (numSeqs/sqrt(2)), (sqrt(3) * numSeqs / 2), distFile + "tempc", cutoff);	
					system(("cat " + distFile + "tempc" + " >> " + distFile).c_str());
					system(("rm " + distFile + "tempc").c_str());				
				}
				else{
					driver(distCalculator, seqDB, (sqrt(3) * numSeqs / 2), numSeqs, distFile + "tempd", cutoff);	
					system(("cat " + distFile + "tempd" + " >> " + distFile).c_str());
					system(("rm " + distFile + "tempd").c_str());				
				}
				wait(NULL);
			}
			wait(NULL);
		}
		wait(NULL);
	
		delete distCalculator;
	
		return 0;

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driver(Dist* distCalculator, SequenceDB* align, int startLine, int endLine, string dFileName, float cutoff){
	try {
		int startTime = time(NULL);
	
		ofstream distFile(dFileName.c_str(), ios::trunc);
		distFile.setf(ios::fixed, ios::showpoint);
		distFile << setprecision(4);
	
		for(int i=startLine;i<endLine;i++){
		
			for(int j=0;j<i;j++){
			
				distCalculator->calcDist(align->get(i), align->get(j));
				double dist = distCalculator->getDist();
				if(dist <= cutoff){
				distFile << align->get(i).getName() << ' ' << align->get(j).getName() << ' ' << dist << endl;
				}
			
			}
			if(i % 100 == 0){
				cout << i << '\t' << time(NULL) - startTime << endl;
			}
		
		}
		cout << endLine-1 << '\t' << time(NULL) - startTime << endl;
	
		return 1;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function driver. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function driver. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	

}

/**************************************************************************************************/

