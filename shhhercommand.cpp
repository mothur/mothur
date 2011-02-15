/*
 *  shhher.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 12/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "shhhercommand.h"

#include "readcolumn.h"
#include "readmatrix.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include "cluster.hpp"
#include "sparsematrix.hpp"
#include <cfloat>

//**********************************************************************************************************************

#define NUMBINS 1000
#define HOMOPS 10
#define MIN_COUNT 0.1
#define MIN_WEIGHT 0.1
#define MIN_TAU 0.0001
#define MIN_ITER 10

//**********************************************************************************************************************

vector<string> ShhherCommand::getValidParameters(){	
	try {
		string Array[] =  {	
			"file", "flow", "lookup", "cutoff", "sigma", "outputdir","inputdir", "processors", "maxiter", "mindelta"	
		};
		
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getValidParameters");
		exit(1);
	}
}

//**********************************************************************************************************************

ShhherCommand::ShhherCommand(){	
	try {
		abort = true; calledHelp = true;
		
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["pn.dist"] = tempOutNames;

	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "ShhherCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<string> ShhherCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"flow"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getRequiredParameters");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<string> ShhherCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getRequiredFiles");
		exit(1);
	}
}

//**********************************************************************************************************************

ShhherCommand::ShhherCommand(string option) {
	try {

#ifdef USE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		MPI_Comm_size(MPI_COMM_WORLD, &ncpus);

		if(pid == 0){
#endif
		
		
		abort = false; calledHelp = false;   
		
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {
				"file", "flow", "lookup", "cutoff", "sigma", "outputdir","inputdir", "processors", "maxiter", "mindelta"	
			};
			
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["pn.dist"] = tempOutNames;
			//			outputTypes["fasta"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
				
				it = parameters.find("lookup");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["lookup"] = inputDir + it->second;		}
				}

				it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			flowFileName = validParameter.validFile(parameters, "flow", true);
			flowFilesFileName = validParameter.validFile(parameters, "file", true);
			if (flowFileName == "not found" && flowFilesFileName == "not found") {
				m->mothurOut("values for either flow or file must be provided for the shhh.seqs command.");
				m->mothurOutEndLine();
				abort = true; 
			}
			else if (flowFileName == "not open" || flowFilesFileName == "not open") { abort = true; }
			
			if(flowFileName != "not found"){	compositeFASTAFileName = "";	}
			else{
				compositeFASTAFileName = flowFilesFileName.substr(0, flowFilesFileName.length()-10) + "pn.fasta";
				ofstream temp;
				m->openOutputFile(compositeFASTAFileName, temp);
				temp.close();
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(flowFileName); //if user entered a file with a path then preserve it	
			}
			
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "lookup", true);
			if (temp == "not found")	{	lookupFileName = "LookUp_Titanium.pat";	}
			else if(temp == "not open")	{	abort = true;			} 
			else						{	lookupFileName = temp;	}
			
			temp = validParameter.validFile(parameters, "processors", false);if (temp == "not found"){	temp = "1";			}
			convert(temp, processors); 

			temp = validParameter.validFile(parameters, "cutoff", false);	if (temp == "not found"){	temp = "0.01";		}
			convert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "mindelta", false);	if (temp == "not found"){	temp = "0.000001";	}
			convert(temp, minDelta); 

			temp = validParameter.validFile(parameters, "maxiter", false);	if (temp == "not found"){	temp = "1000";		}
			convert(temp, maxIters); 

			temp = validParameter.validFile(parameters, "sigma", false);if (temp == "not found")	{	temp = "60";		}
			convert(temp, sigma); 
			
			globaldata = GlobalData::getInstance();
		}
			
#ifdef USE_MPI
		}				
#endif
				
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "ShhherCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

ShhherCommand::~ShhherCommand(){}

//**********************************************************************************************************************

void ShhherCommand::help(){
	try {
		m->mothurOut("The shhher command reads a file containing flowgrams and creates a file of corrected sequences.\n");
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************
#ifdef USE_MPI
int ShhherCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		int tag = 1976;
		MPI_Status status; 

		double begClock = clock();
		unsigned long int begTime = time(NULL);
		
		if(pid == 0){

			for(int i=1;i<ncpus;i++){
				MPI_Send(&abort, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			}
			if(abort == 1){	return 0;	}

			processors = ncpus;
			
			m->mothurOut("\nGetting preliminary data...\n");
			getSingleLookUp();
			getJointLookUp();
			
			vector<string> flowFileVector;
			if(flowFilesFileName != "not found"){
				string fName;

				ifstream flowFilesFile;
				m->openInputFile(flowFilesFileName, flowFilesFile);
				while(flowFilesFile){
					flowFilesFile >> fName;
					flowFileVector.push_back(fName);
					m->gobble(flowFilesFile);
				}
			}
			else{
				flowFileVector.push_back(flowFileName);
			}
			int numFiles = flowFileVector.size();

			for(int i=1;i<ncpus;i++){
				MPI_Send(&numFiles, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			}
			
			for(int i=0;i<numFiles;i++){
				flowFileName = flowFileVector[i];

				
				
				m->mothurOut("\n>>>>>\tProcessing " + flowFileName + " (file " + toString(i+1) + " of " + toString(numFiles) + ")\t<<<<<\n");
				m->mothurOut("Reading flowgrams...\n");
				getFlowData();

				m->mothurOut("Identifying unique flowgrams...\n");
				getUniques();

				m->mothurOut("Calculating distances between flowgrams...\n");
				char fileName[1024];
				strcpy(fileName, flowFileName.c_str());

				for(int i=1;i<ncpus;i++){
					MPI_Send(&fileName[0], 1024, MPI_CHAR, i, tag, MPI_COMM_WORLD);

					MPI_Send(&numSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&numUniques, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&numFlowCells, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&flowDataIntI[0], numSeqs * numFlowCells, MPI_SHORT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&flowDataPrI[0], numSeqs * numFlowCells, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
					MPI_Send(&mapUniqueToSeq[0], numSeqs, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&mapSeqToUnique[0], numSeqs, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&lengths[0], numSeqs, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&jointLookUp[0], NUMBINS * NUMBINS, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
					MPI_Send(&cutoff, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
				}				
							
				string distFileName = flowDistMPI(0, int(sqrt(1.0/float(ncpus)) * numUniques));
				
				int done;
				for(int i=1;i<ncpus;i++){
					MPI_Recv(&done, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
					
					m->appendFiles((distFileName + ".temp." + toString(i)), distFileName);
					remove((distFileName + ".temp." + toString(i)).c_str());
				}

				string namesFileName = createNamesFile();
				
				m->mothurOut("\nClustering flowgrams...\n");
				string listFileName = cluster(distFileName, namesFileName);

				getOTUData(listFileName);
				initPyroCluster();

				for(int i=1;i<ncpus;i++){
					MPI_Send(&numOTUs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&singleLookUp[0], singleLookUp.size(), MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
					MPI_Send(&uniqueFlowgrams[0], numFlowCells * numUniques, MPI_SHORT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&sigma, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
				}
				
				
				double maxDelta = 0;
				int iter = 0;
				
				int numOTUsOnCPU = numOTUs / ncpus;
				int numSeqsOnCPU = numSeqs / ncpus;
				m->mothurOut("\nDenoising flowgrams...\n");
				m->mothurOut("iter\tmaxDelta\tnLL\t\tcycletime\n");
				
				while((maxIters == 0 && maxDelta > minDelta) || iter < MIN_ITER || (maxDelta > minDelta && iter < maxIters)){

					double cycClock = clock();
					unsigned long int cycTime = time(NULL);
					fill();

					int total = singleTau.size();
					for(int i=1;i<ncpus;i++){
						MPI_Send(&total, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&change[0], numOTUs, MPI_SHORT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&centroids[0], numOTUs, MPI_INT, i, tag, MPI_COMM_WORLD);
						
						MPI_Send(&singleTau[0], total, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
						MPI_Send(&seqNumber[0], total, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&seqIndex[0], total, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&nSeqsPerOTU[0], numOTUs, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&cumNumSeqs[0], numOTUs, MPI_INT, i, tag, MPI_COMM_WORLD);
					}
				
					calcCentroidsDriver(0, numOTUsOnCPU);
					
					for(int i=1;i<ncpus;i++){
						int otuStart = i * numOTUs / ncpus;
						int otuStop = (i + 1) * numOTUs / ncpus;
						
						vector<int> tempCentroids(numOTUs, 0);
						vector<short> tempChange(numOTUs, 0);
						
						MPI_Recv(&tempCentroids[0], numOTUs, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
						MPI_Recv(&tempChange[0], numOTUs, MPI_SHORT, i, tag, MPI_COMM_WORLD, &status);
						
						for(int j=otuStart;j<otuStop;j++){
							centroids[j] = tempCentroids[j];
							change[j] = tempChange[j];
						}
					}
									
					maxDelta = getNewWeights();
					double nLL = getLikelihood();
					checkCentroids();
					
					for(int i=1;i<ncpus;i++){
						MPI_Send(&centroids[0], numOTUs, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&weight[0], numOTUs, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
						MPI_Send(&change[0], numOTUs, MPI_SHORT, i, tag, MPI_COMM_WORLD);
					}
					
					calcNewDistancesParent(0, numSeqsOnCPU);

					total = singleTau.size();

					for(int i=1;i<ncpus;i++){
						int childTotal;
						int seqStart = i * numSeqs / ncpus;
						int seqStop = (i + 1) * numSeqs / ncpus;
						
						MPI_Recv(&childTotal, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

						vector<int> childSeqIndex(childTotal, 0);
						vector<double> childSingleTau(childTotal, 0);
						vector<double> childDist(numSeqs * numOTUs, 0);
						vector<int> otuIndex(childTotal, 0);
						
						MPI_Recv(&childSeqIndex[0], childTotal, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
						MPI_Recv(&childSingleTau[0], childTotal, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
						MPI_Recv(&childDist[0], numOTUs * numSeqs, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
						MPI_Recv(&otuIndex[0], childTotal, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

						int oldTotal = total;
						total += childTotal;
						singleTau.resize(total, 0);
						seqIndex.resize(total, 0);
						seqNumber.resize(total, 0);
						
						int childIndex = 0;
						
						for(int j=oldTotal;j<total;j++){
							int otuI = otuIndex[childIndex];
							int seqI = childSeqIndex[childIndex];
							
							singleTau[j] = childSingleTau[childIndex];
							
							aaP[otuI][nSeqsPerOTU[otuI]] = j;
							aaI[otuI][nSeqsPerOTU[otuI]] = seqI;
							nSeqsPerOTU[otuI]++;
							childIndex++;
						}
						
						int index = seqStart * numOTUs;
						for(int j=seqStart;j<seqStop;j++){
							for(int k=0;k<numOTUs;k++){
								dist[index] = childDist[index];
								index++;
							}
						}					
					}
					
					iter++;
					
					m->mothurOut(toString(iter) + '\t' + toString(maxDelta) + '\t' + toString(nLL) + '\t' + toString(time(NULL) - cycTime) + '\t' + toString((clock() - cycClock)/(double)CLOCKS_PER_SEC) + '\n');			

					if((maxIters == 0 && maxDelta > minDelta) || iter < MIN_ITER || (maxDelta > minDelta && iter < maxIters)){
						int live = 1;
						for(int i=1;i<ncpus;i++){
							MPI_Send(&live, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
						}
					}
					else{
						int live = 0;
						for(int i=1;i<ncpus;i++){
							MPI_Send(&live, 1, MPI_INT, i, tag, MPI_COMM_WORLD); //send kill command
						}
					}
					
				}	
				
				m->mothurOut("\nFinalizing...\n");
				fill();
				setOTUs();
				vector<int> otuCounts(numOTUs, 0);
				for(int i=0;i<numSeqs;i++)	{	otuCounts[otuData[i]]++;	}
				calcCentroidsDriver(0, numOTUs);
				writeQualities(otuCounts);
				writeSequences(otuCounts);
				writeNames(otuCounts);
				writeClusters(otuCounts);
				writeGroups();
				
				remove(distFileName.c_str());
				remove(namesFileName.c_str());
				remove(listFileName.c_str());
								 
				m->mothurOut("Total time to process " + toString(flowFileName) + ":\t" + toString(time(NULL) - begTime) + '\t' + toString((clock() - begClock)/(double)CLOCKS_PER_SEC) + '\n');			
			}
		}
		else{
			int abort = 1;
			bool live = 1;

			MPI_Recv(&abort, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			if(abort){	return 0;	}

			int numFiles;
			MPI_Recv(&numFiles, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

			for(int i=0;i<numFiles;i++){
				//Now into the pyrodist part
				char fileName[1024];
				MPI_Recv(&fileName, 1024, MPI_CHAR, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&numSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&numUniques, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&numFlowCells, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				
				flowDataIntI.resize(numSeqs * numFlowCells);
				flowDataPrI.resize(numSeqs * numFlowCells);
				mapUniqueToSeq.resize(numSeqs);
				mapSeqToUnique.resize(numSeqs);
				lengths.resize(numSeqs);
				jointLookUp.resize(NUMBINS * NUMBINS);

				MPI_Recv(&flowDataIntI[0], numSeqs * numFlowCells, MPI_SHORT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&flowDataPrI[0], numSeqs * numFlowCells, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&mapUniqueToSeq[0], numSeqs, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&mapSeqToUnique[0], numSeqs, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&lengths[0], numSeqs, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&jointLookUp[0], NUMBINS * NUMBINS, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&cutoff, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);

				flowFileName = string(fileName);
				int flowDistStart = int(sqrt(float(pid)/float(ncpus)) * numUniques);
				int flowDistEnd = int(sqrt(float(pid+1)/float(ncpus)) * numUniques);
				
				string distanceStringChild = flowDistMPI(flowDistStart, flowDistEnd);

				int done = 1;
				MPI_Send(&done, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);

				//Now into the pyronoise part
				MPI_Recv(&numOTUs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				
				singleLookUp.resize(HOMOPS * NUMBINS);
				uniqueFlowgrams.resize(numUniques * numFlowCells);
				weight.resize(numOTUs);
				centroids.resize(numOTUs);
				change.resize(numOTUs);
				dist.assign(numOTUs * numSeqs, 0);
				nSeqsPerOTU.resize(numOTUs);
				cumNumSeqs.resize(numOTUs);

				MPI_Recv(&singleLookUp[0], singleLookUp.size(), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&uniqueFlowgrams[0], uniqueFlowgrams.size(), MPI_SHORT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(&sigma, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
				
				int startOTU = pid * numOTUs / ncpus;
				int endOTU = (pid + 1) * numOTUs / ncpus;

				int startSeq = pid * numSeqs / ncpus;
				int endSeq = (pid + 1) * numSeqs /ncpus;
				
				int total;

				while(live){

					MPI_Recv(&total, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					singleTau.assign(total, 0.0000);
					seqNumber.assign(total, 0);
					seqIndex.assign(total, 0);
					
					MPI_Recv(&change[0], numOTUs, MPI_SHORT, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&centroids[0], numOTUs, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&singleTau[0], total, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&seqNumber[0], total, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&seqIndex[0], total, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&nSeqsPerOTU[0], total, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&cumNumSeqs[0], numOTUs, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					
					calcCentroidsDriver(startOTU, endOTU);

					MPI_Send(&centroids[0], numOTUs, MPI_INT, 0, tag, MPI_COMM_WORLD);
					MPI_Send(&change[0], numOTUs, MPI_SHORT, 0, tag, MPI_COMM_WORLD);

					
					MPI_Recv(&centroids[0], numOTUs, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&weight[0], numOTUs, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(&change[0], numOTUs, MPI_SHORT, 0, tag, MPI_COMM_WORLD, &status);

					vector<int> otuIndex(total, 0);
					calcNewDistancesChildMPI(startSeq, endSeq, otuIndex);
					total = otuIndex.size();
					
					MPI_Send(&total, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
					MPI_Send(&seqIndex[0], total, MPI_INT, 0, tag, MPI_COMM_WORLD);
					MPI_Send(&singleTau[0], total, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
					MPI_Send(&dist[0], numOTUs * numSeqs, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
					MPI_Send(&otuIndex[0], total, MPI_INT, 0, tag, MPI_COMM_WORLD);
				
					MPI_Recv(&live, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				}
			}
		}		

		MPI_Barrier(MPI_COMM_WORLD);
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "execute");
		exit(1);
	}
}

/**************************************************************************************************/

string ShhherCommand::flowDistMPI(int startSeq, int stopSeq){
	try{		
		ostringstream outStream;
		outStream.setf(ios::fixed, ios::floatfield);
		outStream.setf(ios::dec, ios::basefield);
		outStream.setf(ios::showpoint);
		outStream.precision(6);
		
		int begTime = time(NULL);
		double begClock = clock();
		
		for(int i=startSeq;i<stopSeq;i++){
			for(int j=0;j<i;j++){
				float flowDistance = calcPairwiseDist(mapUniqueToSeq[i], mapUniqueToSeq[j]);
				
				if(flowDistance < 1e-6){
					outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << 0.000000 << endl;
				}
				else if(flowDistance <= cutoff){
					outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << flowDistance << endl;
				}
			}
			if(i % 100 == 0){
				m->mothurOut(toString(i) + '\t' + toString(time(NULL) - begTime) + '\t' + toString((clock()-begClock)/CLOCKS_PER_SEC) + '\n');
			}
		}
		
		m->mothurOut(toString(stopSeq) + '\t' + toString(time(NULL) - begTime) + '\t' + toString((clock()-begClock)/CLOCKS_PER_SEC) + '\n');
		
		string fDistFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.dist";
		if(pid != 0){	fDistFileName += ".temp." + toString(pid);	}

		ofstream distFile(fDistFileName.c_str());
		distFile << outStream.str();		
		distFile.close();
		
		return fDistFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "flowDistParentFork");
		exit(1);
	}
}

#else
//**********************************************************************************************************************

int ShhherCommand::execute(){
	try {
		if (abort == true) { return 0; }
		
		getSingleLookUp();
		getJointLookUp();
				
		vector<string> flowFileVector;
		if(flowFilesFileName != "not found"){
			string fName;
			
			ifstream flowFilesFile;
			m->openInputFile(flowFilesFileName, flowFilesFile);
			while(flowFilesFile){
				flowFilesFile >> fName;
				flowFileVector.push_back(fName);
				m->gobble(flowFilesFile);
			}
		}
		else{
			flowFileVector.push_back(flowFileName);
		}
		int numFiles = flowFileVector.size();
		
		
		for(int i=0;i<numFiles;i++){
			flowFileName = flowFileVector[i];

			m->mothurOut("\n>>>>>\tProcessing " + flowFileName + " (file " + toString(i+1) + " of " + toString(numFiles) + ")\t<<<<<\n");
			m->mothurOut("Reading flowgrams...\n");
			getFlowData();
			
			m->mothurOut("Identifying unique flowgrams...\n");
			getUniques();
			
			
			m->mothurOut("Calculating distances between flowgrams...\n");
			string distFileName = createDistFile(processors);
			string namesFileName = createNamesFile();
				
			m->mothurOut("\nClustering flowgrams...\n");
			string listFileName = cluster(distFileName, namesFileName);
			getOTUData(listFileName);
			
			initPyroCluster();
			
			double maxDelta = 0;
			int iter = 0;
			
			double begClock = clock();
			unsigned long int begTime = time(NULL);

			
			cout << numOTUs << endl;
			
			
			m->mothurOut("\nDenoising flowgrams...\n");
			m->mothurOut("iter\tmaxDelta\tnLL\t\tcycletime\n");
			
			while((maxIters == 0 && maxDelta > minDelta) || iter < MIN_ITER || (maxDelta > minDelta && iter < maxIters)){
				
				double cycClock = clock();
				unsigned long int cycTime = time(NULL);
				fill();
				
				calcCentroids();
				
				maxDelta = getNewWeights();
				double nLL = getLikelihood();
				checkCentroids();
				
				calcNewDistances();

				iter++;
				
				m->mothurOut(toString(iter) + '\t' + toString(maxDelta) + '\t' + toString(nLL) + '\t' + toString(time(NULL) - cycTime) + '\t' + toString((clock() - cycClock)/(double)CLOCKS_PER_SEC) + '\n');

			}	
			
			m->mothurOut("\nFinalizing...\n");
			fill();
			setOTUs();
			
			vector<int> otuCounts(numOTUs, 0);
			for(int i=0;i<numSeqs;i++)	{	otuCounts[otuData[i]]++;	}
			
			calcCentroidsDriver(0, numOTUs);
			writeQualities(otuCounts);
			writeSequences(otuCounts);
			writeNames(otuCounts);
			writeClusters(otuCounts);
			writeGroups();
			
			remove(distFileName.c_str());
			remove(namesFileName.c_str());
			remove(listFileName.c_str());
			
			m->mothurOut("Total time to process " + flowFileName + ":\t" + toString(time(NULL) - begTime) + '\t' + toString((clock() - begClock)/(double)CLOCKS_PER_SEC) + '\n');
		}
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "execute");
		exit(1);
	}
}
#endif
/**************************************************************************************************/

void ShhherCommand::getFlowData(){
	try{
		ifstream flowFile;
		m->openInputFile(flowFileName, flowFile);
		
		string seqName;
		seqNameVector.clear();
		lengths.clear();
		flowDataIntI.clear();
		nameMap.clear();
		
		
		int currentNumFlowCells;
		
		float intensity;
		
		flowFile >> numFlowCells;
		int index = 0;//pcluster
		while(!flowFile.eof()){
			flowFile >> seqName >> currentNumFlowCells;
			lengths.push_back(currentNumFlowCells);

			seqNameVector.push_back(seqName);
			nameMap[seqName] = index++;//pcluster

			for(int i=0;i<numFlowCells;i++){
				flowFile >> intensity;
				if(intensity > 9.99)	{	intensity = 9.99;	}
				int intI = int(100 * intensity + 0.0001);
				flowDataIntI.push_back(intI);
			}
			m->gobble(flowFile);
		}
		flowFile.close();
		
		numSeqs = seqNameVector.size();		
		
		for(int i=0;i<numSeqs;i++){
			int iNumFlowCells = i * numFlowCells;
			for(int j=lengths[i];j<numFlowCells;j++){
				flowDataIntI[iNumFlowCells + j] = 0;
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getFlowData");
		exit(1);
	}
}

/**************************************************************************************************/

void ShhherCommand::getSingleLookUp(){
	try{
		//	these are the -log probabilities that a signal corresponds to a particular homopolymer length
		singleLookUp.assign(HOMOPS * NUMBINS, 0);
		
		int index = 0;
		ifstream lookUpFile;
		m->openInputFile(lookupFileName, lookUpFile);
		
		for(int i=0;i<HOMOPS;i++){
			float logFracFreq;
			lookUpFile >> logFracFreq;
			
			for(int j=0;j<NUMBINS;j++)	{
				lookUpFile >> singleLookUp[index];
				index++;
			}
		}	
		lookUpFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getSingleLookUp");
		exit(1);
	}
}

/**************************************************************************************************/

void ShhherCommand::getJointLookUp(){
	try{
		
		//	the most likely joint probability (-log) that two intenities have the same polymer length
		jointLookUp.resize(NUMBINS * NUMBINS, 0);
		
		for(int i=0;i<NUMBINS;i++){
			for(int j=0;j<NUMBINS;j++){		
				
				double minSum = 100000000;
				
				for(int k=0;k<HOMOPS;k++){
					double sum = singleLookUp[k * NUMBINS + i] + singleLookUp[k * NUMBINS + j];
					
					if(sum < minSum)	{	minSum = sum;		}
				}	
				jointLookUp[i * NUMBINS + j] = minSum;
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getJointLookUp");
		exit(1);
	}
}

/**************************************************************************************************/

double ShhherCommand::getProbIntensity(int intIntensity){                          
	try{

		double minNegLogProb = 100000000; 

		
		for(int i=0;i<HOMOPS;i++){//loop signal strength
			float negLogProb = singleLookUp[i * NUMBINS + intIntensity];
			if(negLogProb < minNegLogProb)	{	minNegLogProb = negLogProb; }
		}
		
		return minNegLogProb;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getProbIntensity");
		exit(1);
	}
}

/**************************************************************************************************/

void ShhherCommand::getUniques(){
	try{
		
		
		numUniques = 0;
		uniqueFlowgrams.assign(numFlowCells * numSeqs, -1);
		uniqueCount.assign(numSeqs, 0);							//	anWeights
		uniqueLengths.assign(numSeqs, 0);
		mapSeqToUnique.assign(numSeqs, -1);
		mapUniqueToSeq.assign(numSeqs, -1);
		
		vector<short> uniqueFlowDataIntI(numFlowCells * numSeqs, -1);
		
		for(int i=0;i<numSeqs;i++){
			int index = 0;
			
			vector<short> current(numFlowCells);
			for(int j=0;j<numFlowCells;j++){
				current[j] = short(((flowDataIntI[i * numFlowCells + j] + 50.0)/100.0));
			}
						
			for(int j=0;j<numUniques;j++){
				int offset = j * numFlowCells;
				bool toEnd = 1;
				
				int shorterLength;
				if(lengths[i] < uniqueLengths[j])	{	shorterLength = lengths[i];			}
				else								{	shorterLength = uniqueLengths[j];	}

				for(int k=0;k<shorterLength;k++){
					if(current[k] != uniqueFlowgrams[offset + k]){
						toEnd = 0;
						break;
					}
				}
				
				if(toEnd){
					mapSeqToUnique[i] = j;
					uniqueCount[j]++;
					index = j;
					if(lengths[i] > uniqueLengths[j])	{	uniqueLengths[j] = lengths[i];	}
					break;
				}
				index++;
			}
			
			if(index == numUniques){
				uniqueLengths[numUniques] = lengths[i];
				uniqueCount[numUniques] = 1;
				mapSeqToUnique[i] = numUniques;//anMap
				mapUniqueToSeq[numUniques] = i;//anF
				
				for(int k=0;k<numFlowCells;k++){
					uniqueFlowgrams[numUniques * numFlowCells + k] = current[k];
					uniqueFlowDataIntI[numUniques * numFlowCells + k] = flowDataIntI[i * numFlowCells + k];
				}
				
				numUniques++;
			}
		}
		uniqueFlowDataIntI.resize(numFlowCells * numUniques);
		uniqueLengths.resize(numUniques);	
		
		flowDataPrI.resize(numSeqs * numFlowCells, 0);
		for(int i=0;i<flowDataPrI.size();i++)	{	flowDataPrI[i] = getProbIntensity(flowDataIntI[i]);		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getUniques");
		exit(1);
	}
}

/**************************************************************************************************/

float ShhherCommand::calcPairwiseDist(int seqA, int seqB){
	try{
		int minLength = lengths[mapSeqToUnique[seqA]];
		if(lengths[seqB] < minLength){	minLength = lengths[mapSeqToUnique[seqB]];	}
		
		int ANumFlowCells = seqA * numFlowCells;
		int BNumFlowCells = seqB * numFlowCells;
		
		float dist = 0;
		
		for(int i=0;i<minLength;i++){
			int flowAIntI = flowDataIntI[ANumFlowCells + i];
			float flowAPrI = flowDataPrI[ANumFlowCells + i];
			
			int flowBIntI = flowDataIntI[BNumFlowCells + i];
			float flowBPrI = flowDataPrI[BNumFlowCells + i];
			dist += jointLookUp[flowAIntI * NUMBINS + flowBIntI] - flowAPrI - flowBPrI;
		}
		
		dist /= (float) minLength;
		return dist;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcPairwiseDist");
		exit(1);
	}
}

/**************************************************************************************************/

void ShhherCommand::flowDistParentFork(string distFileName, int startSeq, int stopSeq){
	try{		

		ostringstream outStream;
		outStream.setf(ios::fixed, ios::floatfield);
		outStream.setf(ios::dec, ios::basefield);
		outStream.setf(ios::showpoint);
		outStream.precision(6);
		
		int begTime = time(NULL);
		double begClock = clock();

		for(int i=startSeq;i<stopSeq;i++){
			for(int j=0;j<i;j++){
				float flowDistance = calcPairwiseDist(mapUniqueToSeq[i], mapUniqueToSeq[j]);

				if(flowDistance < 1e-6){
					outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << 0.000000 << endl;
				}
				else if(flowDistance <= cutoff){
					outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << flowDistance << endl;
				}
			}
			if(i % 100 == 0){
				m->mothurOut(toString(i) + "\t" + toString(time(NULL) - begTime));
				m->mothurOut("\t" + toString((clock()-begClock)/CLOCKS_PER_SEC));
				m->mothurOutEndLine();
			}
		}
		m->mothurOut(toString(stopSeq-1) + "\t" + toString(time(NULL) - begTime));
		m->mothurOut("\t" + toString((clock()-begClock)/CLOCKS_PER_SEC));
		m->mothurOutEndLine();
		
		ofstream distFile(distFileName.c_str());
		distFile << outStream.str();		
		distFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "flowDistParentFork");
		exit(1);
	}
}

/**************************************************************************************************/

string ShhherCommand::createDistFile(int processors){
	try{
		string fDistFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.dist";
				
		unsigned long int begTime = time(NULL);
		double begClock = clock();

		vector<int> start;
		vector<int> end;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors == 1)	{	flowDistParentFork(fDistFileName, 0, numUniques);		}
		else{ //you have multiple processors
			
			if (numSeqs < processors){	processors = 1;	}
			
			vector<int> start(processors, 0);
			vector<int> end(processors, 0);
			
			for (int i = 0; i < processors; i++) {
				start[i] = int(sqrt(float(i)/float(processors)) * numUniques);
				end[i] = int(sqrt(float(i+1)/float(processors)) * numUniques);
			}
			
			int process = 1;
			vector<int> processIDs;
			
			//loop through and create all the processes you want
			while (process != processors) {
				int pid = fork();
				
				if (pid > 0) {
					processIDs.push_back(pid);  //create map from line number to pid so you can append files in correct order later
					process++;
				}else if (pid == 0){
					flowDistParentFork(fDistFileName + toString(getpid()) + ".temp", start[process], end[process]);
					exit(0);
				}else { 
					m->mothurOut("[ERROR]: unable to spawn the necessary processes. Error code: " + toString(pid)); m->mothurOutEndLine(); 
					perror(" : ");
					for (int i=0;i<processIDs.size();i++) {  int temp = processIDs[i]; kill (temp, SIGINT); }
					exit(0);
				}
			}
			
			//parent does its part
			flowDistParentFork(fDistFileName, start[0], end[0]);
			
			//force parent to wait until all the processes are done
			for (int i=0;i<processIDs.size();i++) { 
				int temp = processIDs[i];
				wait(&temp);
			}
			
			//append and remove temp files
			for (int i=0;i<processIDs.size();i++) { 
				m->appendFiles((fDistFileName + toString(processIDs[i]) + ".temp"), fDistFileName);
				remove((fDistFileName + toString(processIDs[i]) + ".temp").c_str());
			}
			
		}
		
#else
		flowDistParentFork(fDistFileName, 0, numUniques);
#endif

		m->mothurOutEndLine();
		
		m->mothurOut("Total time: " + toString(time(NULL) - begTime) + '\t' + toString((clock() - begClock)/CLOCKS_PER_SEC) + '\n');
		

		return fDistFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "createDistFile");
		exit(1);
	}
	
}

/**************************************************************************************************/

string ShhherCommand::createNamesFile(){
	try{
		
		vector<string> duplicateNames(numUniques, "");
		for(int i=0;i<numSeqs;i++){
			duplicateNames[mapSeqToUnique[i]] += seqNameVector[i] + ',';
		}
		
		string nameFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.names";
		
		ofstream nameFile;
		m->openOutputFile(nameFileName, nameFile);
		
		for(int i=0;i<numUniques;i++){
//			nameFile << seqNameVector[mapUniqueToSeq[i]] << '\t' << duplicateNames[i].substr(0, duplicateNames[i].find_last_of(',')) << endl;
			nameFile << mapUniqueToSeq[i] << '\t' << duplicateNames[i].substr(0, duplicateNames[i].find_last_of(',')) << endl;
		}
		
		nameFile.close();
		return  nameFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "createNamesFile");
		exit(1);
	}
}

//**********************************************************************************************************************

string ShhherCommand::cluster(string distFileName, string namesFileName){
	try {
		
		
		globaldata->setNameFile(namesFileName);
		globaldata->setColumnFile(distFileName);
		globaldata->setFormat("column");
		
		ReadMatrix* read = new ReadColumnMatrix(distFileName); 	
		read->setCutoff(cutoff);
		
		NameAssignment* clusterNameMap = new NameAssignment(namesFileName);
		clusterNameMap->readMap();
		read->read(clusterNameMap);
		
		ListVector* list = read->getListVector();
		SparseMatrix* matrix = read->getMatrix();
		
		delete read; 
		delete clusterNameMap; 
				
		RAbundVector* rabund = new RAbundVector(list->getRAbundVector());
		
		Cluster* cluster = new CompleteLinkage(rabund, list, matrix, cutoff, "furthest"); 
		string tag = cluster->getTag();
		
		double clusterCutoff = cutoff;
		while (matrix->getSmallDist() <= clusterCutoff && matrix->getNNodes() > 0){
			cluster->update(clusterCutoff);
		}
		
		list->setLabel(toString(cutoff));
		
		string listFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.list";
		ofstream listFile;
		m->openOutputFile(listFileName, listFile);
		list->print(listFile);
		listFile.close();
		
		delete matrix;	delete cluster;	delete rabund; delete list;
	
		return listFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "cluster");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::getOTUData(string listFileName){
	try {

		ifstream listFile;
		m->openInputFile(listFileName, listFile);
		string label;
		
		listFile >> label >> numOTUs;

		otuData.assign(numSeqs, 0);
		cumNumSeqs.assign(numOTUs, 0);
		nSeqsPerOTU.assign(numOTUs, 0);
		aaP.clear();aaP.resize(numOTUs);
		
		seqNumber.clear();
		aaI.clear();
		seqIndex.clear();
		
		string singleOTU = "";
		
		for(int i=0;i<numOTUs;i++){

			listFile >> singleOTU;
			
			istringstream otuString(singleOTU);

			while(otuString){
				
				string seqName = "";
				
				for(int j=0;j<singleOTU.length();j++){
					char letter = otuString.get();
					
					if(letter != ','){
						seqName += letter;
					}
					else{
						map<string,int>::iterator nmIt = nameMap.find(seqName);
						int index = nmIt->second;
						
						nameMap.erase(nmIt);
						
						otuData[index] = i;
						nSeqsPerOTU[i]++;
						aaP[i].push_back(index);
						seqName = "";
					}
				}
				
				map<string,int>::iterator nmIt = nameMap.find(seqName);

				int index = nmIt->second;
				nameMap.erase(nmIt);

				otuData[index] = i;
				nSeqsPerOTU[i]++;
				aaP[i].push_back(index);	
				
				otuString.get();
			}
			
			sort(aaP[i].begin(), aaP[i].end());
			for(int j=0;j<nSeqsPerOTU[i];j++){
				seqNumber.push_back(aaP[i][j]);
			}
			for(int j=nSeqsPerOTU[i];j<numSeqs;j++){
				aaP[i].push_back(0);
			}
			
			
		}
		
		for(int i=1;i<numOTUs;i++){
			cumNumSeqs[i] = cumNumSeqs[i-1] + nSeqsPerOTU[i-1];
		}
		aaI = aaP;
		seqIndex = seqNumber;
		
		listFile.close();	
		
		for(int i=0;i<seqNumber.size();i++){
			cout << seqNumber[i] << ' ';
		}
		cout << endl;
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getOTUData");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::initPyroCluster(){                          
	try{
		dist.assign(numSeqs * numOTUs, 0);
		change.assign(numOTUs, 1);
		centroids.assign(numOTUs, -1);
		weight.assign(numOTUs, 0);
		singleTau.assign(numSeqs, 1.0);
		
		nSeqsBreaks.assign(processors+1, 0);
		nOTUsBreaks.assign(processors+1, 0);

		cout << numSeqs << '\t' << numOTUs << '\t' << processors << endl;
		
		nSeqsBreaks[0] = 0;
		for(int i=0;i<processors;i++){
			nSeqsBreaks[i+1] = nSeqsBreaks[i] + (int)((double) numSeqs / (double) processors);
			nOTUsBreaks[i+1] = nOTUsBreaks[i] + (int)((double) numOTUs / (double) processors);
		}
		nSeqsBreaks[processors] = numSeqs;
		nOTUsBreaks[processors] = numOTUs;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "initPyroCluster");
		exit(1);	
	}
}

/**************************************************************************************************/

void ShhherCommand::fill(){
	try {
		int index = 0;
		for(int i=0;i<numOTUs;i++){
			cumNumSeqs[i] = index;
			for(int j=0;j<nSeqsPerOTU[i];j++){
				seqNumber[index] = aaP[i][j];
				seqIndex[index] = aaI[i][j];
				
				index++;
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "fill");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::calcCentroids(){                          
	try{
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)

		if(processors == 1)	{
			calcCentroidsDriver(0, numOTUs);		
		}
		else{ //you have multiple processors
			if (numOTUs < processors){	processors = 1;	}
			
			int process = 1;
			vector<int> processIDs;
			
			//loop through and create all the processes you want
			while (process != processors) {
				int pid = vfork();
				
				if (pid > 0) {
					processIDs.push_back(pid);  //create map from line number to pid so you can append files in correct order later
					process++;
				}else if (pid == 0){
					calcCentroidsDriver(nOTUsBreaks[process], nOTUsBreaks[process+1]);
					exit(0);
				}else { 
					m->mothurOut("[ERROR]: unable to spawn the necessary processes. Error code: " + toString(pid)); m->mothurOutEndLine(); 
					perror(" : ");
					for (int i=0;i<processIDs.size();i++) {  int temp = processIDs[i]; kill (temp, SIGINT); }
					exit(0);
				}
			}
			
			//parent does its part
			calcCentroidsDriver(nOTUsBreaks[0], nOTUsBreaks[1]);

			//force parent to wait until all the processes are done
			for (int i=0;i<processIDs.size();i++) { 
				int temp = processIDs[i];
				wait(&temp);
			}
		}
		
#else
		calcCentroidsDriver(0, numOTUs);
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcCentroidsDriver");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::calcCentroidsDriver(int start, int finish){                          
	
	//this function gets the most likely homopolymer length at a flow position for a group of sequences
	//within an otu
	
	try{
		
//		for(int i=0;i<seqNumber.size();i++){
//			cout << seqNumber[i] << ' ';
//		}cout << endl;
		
		for(int i=start;i<finish;i++){
			
			double count = 0;
			int position = 0;
			int minFlowGram = 100000000;
			double minFlowValue = 1e8;
			change[i] = 0; //FALSE
			
			for(int j=0;j<nSeqsPerOTU[i];j++){
				count += singleTau[seqNumber[cumNumSeqs[i] + j]];
			}

			if(nSeqsPerOTU[i] > 0 && count > MIN_COUNT){
				vector<double> adF(nSeqsPerOTU[i]);
				vector<int> anL(nSeqsPerOTU[i]);
				
				for(int j=0;j<nSeqsPerOTU[i];j++){
					int index = cumNumSeqs[i] + j;
					int nI = seqIndex[index];
					int nIU = mapSeqToUnique[nI];
					
					int k;
					for(k=0;k<position;k++){
						if(nIU == anL[k]){
							break;
						}
					}
					if(k == position){
						anL[position] = nIU;
						adF[position] = 0.0000;
						position++;
					}						
				}
				
				for(int j=0;j<nSeqsPerOTU[i];j++){
					int index = cumNumSeqs[i] + j;
					int nI = seqIndex[index];
					
					double tauValue = singleTau[seqNumber[index]];
					
					for(int k=0;k<position;k++){
						double dist = getDistToCentroid(anL[k], nI, lengths[nI]);
						adF[k] += dist * tauValue;
					}
				}
				
				for(int j=0;j<position;j++){
					if(adF[j] < minFlowValue){
						minFlowGram = j;
						minFlowValue = adF[j];
					}
				}
				
				if(centroids[i] != anL[minFlowGram]){
					change[i] = 1;
					centroids[i] = anL[minFlowGram];
				}
			}
			else if(centroids[i] != -1){
				change[i] = 1;
				centroids[i] = -1;			
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcCentroidsDriver");
		exit(1);	
	}		
}

/**************************************************************************************************/

double ShhherCommand::getDistToCentroid(int cent, int flow, int length){
	try{
		
		int flowAValue = cent * numFlowCells;
		int flowBValue = flow * numFlowCells;
		
		double dist = 0;
		
		for(int i=0;i<length;i++){
			dist += singleLookUp[uniqueFlowgrams[flowAValue] * NUMBINS + flowDataIntI[flowBValue]];
			flowAValue++;
			flowBValue++;
		}
		
		return dist / (double)length;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getDistToCentroid");
		exit(1);	
	}		
}

/**************************************************************************************************/

double ShhherCommand::getNewWeights(){
	try{
		
		double maxChange = 0;
		
		for(int i=0;i<numOTUs;i++){
			
			double difference = weight[i];
			weight[i] = 0;
			
			for(int j=0;j<nSeqsPerOTU[i];j++){
				int index = cumNumSeqs[i] + j;
				double tauValue = singleTau[seqNumber[index]];
				weight[i] += tauValue;
			}
			
			difference = fabs(weight[i] - difference);
			if(difference > maxChange){	maxChange = difference;	}
		}
		return maxChange;
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getNewWeights");
		exit(1);	
	}		
}

/**************************************************************************************************/

double ShhherCommand::getLikelihood(){
	
	try{
		
		vector<long double> P(numSeqs, 0);
		int effNumOTUs = 0;
		
		for(int i=0;i<numOTUs;i++){
			if(weight[i] > MIN_WEIGHT){
				effNumOTUs++;
			}
		}
		
		string hold;
		for(int i=0;i<numOTUs;i++){
			for(int j=0;j<nSeqsPerOTU[i];j++){
				int index = cumNumSeqs[i] + j;
				int nI = seqIndex[index];
				double singleDist = dist[seqNumber[index]];
				
				P[nI] += weight[i] * exp(-singleDist * sigma);
			}
		}
		double nLL = 0.00;
		for(int i=0;i<numSeqs;i++){
			if(P[i] == 0){	P[i] = DBL_EPSILON;	}

			nLL += -log(P[i]);
		}
		
		nLL = nLL -(double)numSeqs * log(sigma);

		return nLL; 
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "getNewWeights");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::checkCentroids(){
	try{
		vector<int> unique(numOTUs, 1);
		
		for(int i=0;i<numOTUs;i++){
			if(centroids[i] == -1 || weight[i] < MIN_WEIGHT){
				unique[i] = -1;
			}
		}
		
		for(int i=0;i<numOTUs;i++){
			if(unique[i] == 1){
				for(int j=i+1;j<numOTUs;j++){
					if(unique[j] == 1){
						
						if(centroids[j] == centroids[i]){
							unique[j] = 0;
							centroids[j] = -1;
							
							weight[i] += weight[j];
							weight[j] = 0.0;
						}
					}
				}
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "checkCentroids");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::calcNewDistances(){                          
	try{
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		
		if(processors == 1)	{
			calcNewDistancesParent(0, numSeqs);		
		}
		else{ //you have multiple processors
			if (numSeqs < processors){	processors = 1;	}
			
			vector<vector<int> > child_otuIndex(processors);
			vector<vector<int> > child_seqIndex(processors);
			vector<vector<double> > child_singleTau(processors);			
			vector<int> totals(processors);
			
			int process = 1;
			vector<int> processIDs;

			//loop through and create all the processes you want
			while (process != processors) {
				int pid = vfork();
				
				if (pid > 0) {
					processIDs.push_back(pid);  //create map from line number to pid so you can append files in correct order later
					process++;
				}else if (pid == 0){
					calcNewDistancesChild(nSeqsBreaks[process], nSeqsBreaks[process+1], child_otuIndex[process], child_seqIndex[process], child_singleTau[process]);
					totals[process] = child_otuIndex[process].size();

					exit(0);
				}else { 
					m->mothurOut("[ERROR]: unable to spawn the necessary processes. Error code: " + toString(pid)); m->mothurOutEndLine(); 
					perror(" : ");
					for (int i=0;i<processIDs.size();i++) {  int temp = processIDs[i]; kill (temp, SIGINT); }
					exit(0);
				}
			}
			
			//parent does its part
			calcNewDistancesParent(nSeqsBreaks[0], nSeqsBreaks[1]);
			int total = seqIndex.size();

			//force parent to wait until all the processes are done
			for (int i=0;i<processIDs.size();i++) { 
				int temp = processIDs[i];
				wait(&temp);
			}

			for(int i=1;i<processors;i++){
				int oldTotal = total;
				total += totals[i];

				singleTau.resize(total, 0);
				seqIndex.resize(total, 0);
				seqNumber.resize(total, 0);
				
				int childIndex = 0;
				
				for(int j=oldTotal;j<total;j++){
					int otuI = child_otuIndex[i][childIndex];
					int seqI = child_seqIndex[i][childIndex];

					singleTau[j] = child_singleTau[i][childIndex];
					aaP[otuI][nSeqsPerOTU[otuI]] = j;
					aaI[otuI][nSeqsPerOTU[otuI]] = seqI;
					nSeqsPerOTU[otuI]++;

					childIndex++;
				}
			}
		}
#else
		calcNewDistancesParent(0, numSeqs);		
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcNewDistances");
		exit(1);	
	}		
}

/**************************************************************************************************/
#ifdef USE_MPI
void ShhherCommand::calcNewDistancesChildMPI(int startSeq, int stopSeq, vector<int>& otuIndex){
	
	try{
		vector<double> newTau(numOTUs,0);
		vector<double> norms(numSeqs, 0);
		otuIndex.clear();
		seqIndex.clear();
		singleTau.clear();
		
		
		
		for(int i=startSeq;i<stopSeq;i++){
			double offset = 1e8;
			int indexOffset = i * numOTUs;
			
			for(int j=0;j<numOTUs;j++){
				
				if(weight[j] > MIN_WEIGHT && change[j] == 1){
					dist[indexOffset + j] = getDistToCentroid(centroids[j], i, lengths[i]);
				}
				if(weight[j] > MIN_WEIGHT && dist[indexOffset + j] < offset){
					offset = dist[indexOffset + j];
				}
			}
			
			for(int j=0;j<numOTUs;j++){
				if(weight[j] > MIN_WEIGHT){
					newTau[j] = exp(sigma * (-dist[indexOffset + j] + offset)) * weight[j];
					norms[i] += newTau[j];
				}
				else{
					newTau[j] = 0.0;
				}
			}
			
			for(int j=0;j<numOTUs;j++){

				newTau[j] /= norms[i];
				
				if(newTau[j] > MIN_TAU){
					otuIndex.push_back(j);
					seqIndex.push_back(i);
					singleTau.push_back(newTau[j]);
				}
			}
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcNewDistancesChildMPI");
		exit(1);	
	}		
}
#endif
/**************************************************************************************************/

void ShhherCommand::calcNewDistancesChild(int startSeq, int stopSeq, vector<int>& child_otuIndex, vector<int>& child_seqIndex, vector<double>& child_singleTau){
	
	try{
		vector<double> newTau(numOTUs,0);
		vector<double> norms(numSeqs, 0);
		child_otuIndex.resize(0);
		child_seqIndex.resize(0);
		child_singleTau.resize(0);
		
		for(int i=startSeq;i<stopSeq;i++){
			double offset = 1e8;
			int indexOffset = i * numOTUs;
			
			
			for(int j=0;j<numOTUs;j++){
				if(weight[j] > MIN_WEIGHT && change[j] == 1){
					dist[indexOffset + j] = getDistToCentroid(centroids[j], i, lengths[i]);
				}
				
				if(weight[j] > MIN_WEIGHT && dist[indexOffset + j] < offset){
					offset = dist[indexOffset + j];
				}
			}
			
			for(int j=0;j<numOTUs;j++){
				if(weight[j] > MIN_WEIGHT){
					newTau[j] = exp(sigma * (-dist[indexOffset + j] + offset)) * weight[j];
					norms[i] += newTau[j];
				}
				else{
					newTau[j] = 0.0;
				}
			}
			
			for(int j=0;j<numOTUs;j++){
				newTau[j] /= norms[i];
				
				if(newTau[j] > MIN_TAU){
					child_otuIndex.push_back(j);
					child_seqIndex.push_back(i);
					child_singleTau.push_back(newTau[j]);
				}
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcNewDistancesChild");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::calcNewDistancesParent(int startSeq, int stopSeq){
	
	try{
		
		int total = 0;
		vector<double> newTau(numOTUs,0);
		vector<double> norms(numSeqs, 0);
		nSeqsPerOTU.assign(numOTUs, 0);
		
		for(int i=startSeq;i<stopSeq;i++){
			int indexOffset = i * numOTUs;
			
			double offset = 1e8;
			
			for(int j=0;j<numOTUs;j++){
				if(weight[j] > MIN_WEIGHT && change[j] == 1){
					dist[indexOffset + j] = getDistToCentroid(centroids[j], i, lengths[i]);
				}
				
				if(weight[j] > MIN_WEIGHT && dist[indexOffset + j] < offset){
					offset = dist[indexOffset + j];
				}
			}
			
			for(int j=0;j<numOTUs;j++){
				if(weight[j] > MIN_WEIGHT){
					newTau[j] = exp(sigma * (-dist[indexOffset + j] + offset)) * weight[j];
					norms[i] += newTau[j];
				}
				else{
					newTau[j] = 0.0;
				}
			}
			
			for(int j=0;j<numOTUs;j++){
				newTau[j] /= norms[i];
			}
			
			for(int j=0;j<numOTUs;j++){
				if(newTau[j] > MIN_TAU){
					
					int oldTotal = total;
					
					total++;
					
					singleTau.resize(total, 0);
					seqNumber.resize(total, 0);
					seqIndex.resize(total, 0);
					
					singleTau[oldTotal] = newTau[j];
					
					aaP[j][nSeqsPerOTU[j]] = oldTotal;
					aaI[j][nSeqsPerOTU[j]] = i;
					nSeqsPerOTU[j]++;
				}
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcNewDistancesParent");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::setOTUs(){
	
	try {
		vector<double> bigTauMatrix(numOTUs * numSeqs, 0.0000);
		
		for(int i=0;i<numOTUs;i++){
			for(int j=0;j<nSeqsPerOTU[i];j++){
				int index = cumNumSeqs[i] + j;
				double tauValue = singleTau[seqNumber[index]];
				int sIndex = seqIndex[index];
				bigTauMatrix[sIndex * numOTUs + i] = tauValue;				
			}
		}
		
		for(int i=0;i<numSeqs;i++){
			double maxTau = -1.0000;
			int maxOTU = -1;
			for(int j=0;j<numOTUs;j++){
				if(bigTauMatrix[i * numOTUs + j] > maxTau){
					maxTau = bigTauMatrix[i * numOTUs + j];
					maxOTU = j;
				}
			}
			
			otuData[i] = maxOTU;
		}
		
		nSeqsPerOTU.assign(numOTUs, 0);		
		
		for(int i=0;i<numSeqs;i++){
			int index = otuData[i];
			
			singleTau[i] = 1.0000;
			dist[i] = 0.0000;
			
			aaP[index][nSeqsPerOTU[index]] = i;
			aaI[index][nSeqsPerOTU[index]] = i;
			
			nSeqsPerOTU[index]++;
		}
		fill();	
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "calcNewDistances");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeQualities(vector<int> otuCounts){
	
	try {
		string qualityFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.qual";

		ofstream qualityFile;
		m->openOutputFile(qualityFileName, qualityFile);

		qualityFile.setf(ios::fixed, ios::floatfield);
		qualityFile.setf(ios::showpoint);
		qualityFile << setprecision(6);
		
		vector<vector<int> > qualities(numOTUs);
		vector<double> pr(HOMOPS, 0);
		
		
		for(int i=0;i<numOTUs;i++){
			int index = 0;
			int base = 0;
			
			if(nSeqsPerOTU[i] > 0){
				qualities[i].assign(1024, -1);
				
				while(index < numFlowCells){
					double maxPrValue = 1e8;
					short maxPrIndex = -1;
					double count = 0.0000;
					
					pr.assign(HOMOPS, 0);
					
					for(int j=0;j<nSeqsPerOTU[i];j++){
						int lIndex = cumNumSeqs[i] + j;
						double tauValue = singleTau[seqNumber[lIndex]];
						int sequenceIndex = aaI[i][j];
						short intensity = flowDataIntI[sequenceIndex * numFlowCells + index];
						
						count += tauValue;
						
						for(int s=0;s<HOMOPS;s++){
							pr[s] += tauValue * singleLookUp[s * NUMBINS + intensity];
						}
					}
					
					maxPrIndex = uniqueFlowgrams[centroids[i] * numFlowCells + index];
					maxPrValue = pr[maxPrIndex];
					
					if(count > MIN_COUNT){
						double U = 0.0000;
						double norm = 0.0000;
						
						for(int s=0;s<HOMOPS;s++){
							norm += exp(-(pr[s] - maxPrValue));
						}
						
						for(int s=1;s<=maxPrIndex;s++){
							int value = 0;
							double temp = 0.0000;
							
							U += exp(-(pr[s-1]-maxPrValue))/norm;
							
							if(U>0.00){
								temp = log10(U);
							}
							else{
								temp = -10.1;
							}
							temp = floor(-10 * temp);
							value = (int)floor(temp);
							if(value > 100){	value = 100;	}
							
							qualities[i][base] = (int)value;
							base++;
						}
					}
					
					index++;
				}
			}
			
			
			if(otuCounts[i] > 0){
				qualityFile << '>' << seqNameVector[mapUniqueToSeq[i]] << endl;
				
				int j=4;	//need to get past the first four bases
				while(qualities[i][j] != -1){
					qualityFile << qualities[i][j] << ' ';
					j++;
				}
				qualityFile << endl;
			}
		}
		qualityFile.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeQualities");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeSequences(vector<int> otuCounts){
	try {
		string bases = "TACG";
		
		string fastaFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.fasta";
		ofstream fastaFile;
		m->openOutputFile(fastaFileName, fastaFile);
		
		vector<string> names(numOTUs, "");
		
		for(int i=0;i<numOTUs;i++){
			int index = centroids[i];
			
			if(otuCounts[i] > 0){
				fastaFile << '>' << seqNameVector[aaI[i][0]] << endl;
				
				for(int j=8;j<numFlowCells;j++){
					
					char base = bases[j % 4];
					for(int k=0;k<uniqueFlowgrams[index * numFlowCells + j];k++){
						fastaFile << base;
					}
				}
				fastaFile << endl;
			}
		}
		fastaFile.close();
		
		if(compositeFASTAFileName != ""){
			m->appendFiles(fastaFileName, compositeFASTAFileName);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeSequences");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeNames(vector<int> otuCounts){
	try {
		string nameFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.final.names";
		ofstream nameFile;
		m->openOutputFile(nameFileName, nameFile);
		
		for(int i=0;i<numOTUs;i++){
			if(otuCounts[i] > 0){
				nameFile << seqNameVector[aaI[i][0]] << '\t' << seqNameVector[aaI[i][0]];
				
				for(int j=1;j<nSeqsPerOTU[i];j++){
					nameFile << ',' << seqNameVector[aaI[i][j]];
				}
				
				nameFile << endl;
			}
		}
		nameFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeNames");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeGroups(){
	try {
		string fileRoot = flowFileName.substr(0,flowFileName.find_last_of('.'));
		string groupFileName = fileRoot + ".pn.groups";
		ofstream groupFile;
		m->openOutputFile(groupFileName, groupFile);
		
		for(int i=0;i<numSeqs;i++){
			groupFile << seqNameVector[i] << '\t' << fileRoot << endl;
		}
		groupFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeGroups");
		exit(1);	
	}		
}

/**************************************************************************************************/

void ShhherCommand::writeClusters(vector<int> otuCounts){
	try {
		string otuCountsFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".pn.counts";
		ofstream otuCountsFile;
		m->openOutputFile(otuCountsFileName, otuCountsFile);
		
		string bases = "TACG";
		
		for(int i=0;i<numOTUs;i++){
			//output the translated version of the centroid sequence for the otu
			if(otuCounts[i] > 0){
				int index = centroids[i];
				
				otuCountsFile << "ideal\t";
				for(int j=8;j<numFlowCells;j++){
					char base = bases[j % 4];
					for(int s=0;s<uniqueFlowgrams[index * numFlowCells + j];s++){
						otuCountsFile << base;
					}
				}
				otuCountsFile << endl;
				
				for(int j=0;j<nSeqsPerOTU[i];j++){
					int sequence = aaI[i][j];
					otuCountsFile << seqNameVector[sequence] << '\t';
					
					for(int k=8;k<lengths[sequence];k++){
						char base = bases[k % 4];
						int freq = int(0.01 * (double)flowDataIntI[sequence * numFlowCells + k] + 0.5);
						
						for(int s=0;s<freq;s++){
							otuCountsFile << base;
						}
					}
					otuCountsFile << endl;
				}
				otuCountsFile << endl;
			}
		}
		otuCountsFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "writeClusters");
		exit(1);	
	}		
}

//**********************************************************************************************************************
