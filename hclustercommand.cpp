/*
 *  hclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/13/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "hclustercommand.h"

//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
HClusterCommand::HClusterCommand(string option){
	try{
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"cutoff","precision","method","showabund","timing","phylip","column","name","sorted"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {
					abort = true;
				}
			}
			
			globaldata->newRead();
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  distfile = phylipfile;  format = "phylip"; 	}
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  distfile = columnfile; format = "column";	}
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			
			if ((phylipfile == "") && (columnfile == "")) { mothurOut("When executing a cluster command you must enter a phylip or a column."); mothurOutEndLine(); abort = true; }
			else if ((phylipfile != "") && (columnfile != "")) { mothurOut("When executing a cluster command you must enter ONLY ONE of the following: phylip or column."); mothurOutEndLine(); abort = true; }
		
			if (columnfile != "") {
				if (namefile == "") {  cout << "You need to provide a namefile if you are going to use the column format." << endl; abort = true; }
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			//get user cutoff and precision or use defaults
			string temp;
			temp = validParameter.validFile(parameters, "precision", false);
			if (temp == "not found") { temp = "100"; }
			//saves precision legnth for formatting below
			length = temp.length();
			convert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "nearest"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else { mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average."); mothurOutEndLine(); abort = true; }

			showabund = validParameter.validFile(parameters, "showabund", false);
			if (showabund == "not found") { showabund = "T"; }
			
			sort = validParameter.validFile(parameters, "sorted", false);
			if (sort == "not found") { sort = "F"; }
			sorted = isTrue(sort);

			timing = validParameter.validFile(parameters, "timing", false);
			if (timing == "not found") { timing = "F"; }
			
				
			if (abort == false) {
											
				fileroot = getRootName(distfile);
				
				tag = "fn";  //until we figure out average and nearest methods
			
				openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
				openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
				openOutputFile(fileroot+ tag + ".list",		listFile);
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "HClusterCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void HClusterCommand::help(){
	try {
		mothurOut("The cluster command can only be executed after a successful read.dist command.\n");
		mothurOut("The cluster command parameter options are method, cuttoff, precision, showabund and timing. No parameters are required.\n");
		mothurOut("The cluster command should be in the following format: \n");
		mothurOut("cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n");
		mothurOut("The acceptable cluster methods are furthest, nearest and average.  If no method is provided then furthest is assumed.\n\n");	
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

HClusterCommand::~HClusterCommand(){}

//**********************************************************************************************************************

int HClusterCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		if(namefile != ""){	
			globaldata->nameMap = new NameAssignment(namefile);
			globaldata->nameMap->readMap();
		}else{
			globaldata->nameMap = NULL;
		}
		
		time_t estart = time(NULL);
		
		if (!sorted) {
			read = new ReadCluster(distfile, cutoff); 	
			read->setFormat(format);
			read->read(globaldata->nameMap);
			distfile = read->getOutputFile();
		
			list = read->getListVector();
			delete read;
		}else {
			list = new ListVector(globaldata->nameMap->getListVector());
		}
	
		mothurOut("It took " + toString(time(NULL) - estart) + " seconds to sort. "); mothurOutEndLine();
		estart = time(NULL);
		
		//list vector made by read contains all sequence names
		if(list != NULL){
			rabund = new RAbundVector(list->getRAbundVector());
		}else{
			mothurOut("Error: no list vector!"); mothurOutEndLine(); return 0;
		}
		
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;
		
		print_start = true;
		start = time(NULL);
		
//cout << "here" << endl;	
		ifstream in;
		openInputFile(distfile, in);
		string firstName, secondName;
		float distance;
		
		cluster = new HCluster(rabund, list);
		bool clusteredSomething;
		vector<seqDist> seqs; seqs.resize(1); // to start loop
		exitedBreak = false;  //lets you know if there is a distance stored in next
		
		while (seqs.size() != 0){
		
			seqs = getSeqs(in);
			random_shuffle(seqs.begin(), seqs.end());
			
			if (seqs.size() == 0) { break; } //there are no more distances
		
			for (int i = 0; i < seqs.size(); i++) {  //-1 means skip me

				if (print_start && isTrue(timing)) {
					mothurOut("Clustering (" + tag + ") dist " + toString(distance) + "/" 
							  + toString(roundDist(distance, precision)) 
							  + "\t(precision: " + toString(precision) + ")");
					cout.flush();
					print_start = false;
				}
				
	///cout << "before cluster update" << endl;
				if (seqs[i].seq1 != seqs[i].seq2) {
					clusteredSomething = cluster->update(seqs[i].seq1, seqs[i].seq2, seqs[i].dist);
					
					float rndDist = roundDist(seqs[i].dist, precision);
					//cout << "after cluster update clusterSomething = " << clusteredSomething << " rndDist = " << rndDist << " rndPreviousDist = " << rndPreviousDist << endl;			
					
					
					if((previousDist <= 0.0000) && (seqs[i].dist != previousDist)){
						printData("unique");
					}
					else if((rndDist != rndPreviousDist)){
						printData(toString(rndPreviousDist,  length-1));
					}
					
					previousDist = seqs[i].dist;
					rndPreviousDist = rndDist;
					oldRAbund = *rabund;
					oldList = *list;
				}
			}
		}
		
		in.close();

		if (print_start && isTrue(timing)) {
			//mothurOut("Clustering (" + tag + ") for distance " + toString(previousDist) + "/" + toString(rndPreviousDist) 
					 //+ "\t(precision: " + toString(precision) + ", Nodes: " + toString(matrix->getNNodes()) + ")");
			cout.flush();
	 		print_start = false;
		}
	
		if(previousDist <= 0.0000){
			printData("unique");
		}
		else if(rndPreviousDist<cutoff){
			printData(toString(rndPreviousDist, length-1));
		}
		
		//delete globaldata's copy of the sparsematrix and listvector to free up memory
		delete globaldata->gListVector;	 globaldata->gListVector = NULL;
		
		//saves .list file so you can do the collect, rarefaction and summary commands without doing a read.list
		if (globaldata->getFormat() == "phylip") { globaldata->setPhylipFile(""); }
		else if (globaldata->getFormat() == "column") { globaldata->setColumnFile(""); }
		
		globaldata->setListFile(fileroot+ tag + ".list");
		globaldata->setNameFile("");
		globaldata->setFormat("list");
		
		sabundFile.close();
		rabundFile.close();
		listFile.close();
		
		delete cluster;
		//if (isTrue(timing)) {
			mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster. "); mothurOutEndLine();
		//}
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void HClusterCommand::printData(string label){
	try {
		if (isTrue(timing)) {
			mothurOut("\tTime: " + toString(time(NULL) - start) + "\tsecs for " + toString(oldRAbund.getNumBins()) 
		     + "\tclusters. Updates: " + toString(loops)); mothurOutEndLine();
		}
		print_start = true;
		loops = 0;
		start = time(NULL);

		oldRAbund.setLabel(label);
		if (isTrue(showabund)) {
			oldRAbund.getSAbundVector().print(cout);
		}
		oldRAbund.print(rabundFile);
		oldRAbund.getSAbundVector().print(sabundFile);
	
		oldList.setLabel(label);
		oldList.print(listFile);
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "printData");
		exit(1);
	}


}
//**********************************************************************************************************************

vector<seqDist> HClusterCommand::getSeqs(ifstream& filehandle){
	try {
		string firstName, secondName;
		float distance, prevDistance;
		vector<seqDist> sameSeqs;
		prevDistance = -1;
		
		//if you are not at the beginning of the file
		if (exitedBreak) { 
			sameSeqs.push_back(next);
			prevDistance = next.dist;
			exitedBreak = false;
		}
	
		//get entry
		while (filehandle) {
			
			filehandle >> firstName >> secondName >> distance;  
//cout << firstName << '\t' << secondName << '\t' << distance << endl;
			gobble(filehandle);
			
			//save first one
			if (prevDistance == -1) { prevDistance = distance; }
			
			map<string,int>::iterator itA = globaldata->nameMap->find(firstName);
			map<string,int>::iterator itB = globaldata->nameMap->find(secondName);
			
			if(itA == globaldata->nameMap->end()){
				cerr << "AAError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
			}
			if(itB == globaldata->nameMap->end()){
				cerr << "ABError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
			}
			
			//using cutoff
			if (distance > cutoff) { break; }
			
			if (distance != -1) { //-1 means skip me
				
				//are the distances the same
				if (distance == prevDistance) { //save in vector
					seqDist temp;
					temp.seq1 = itA->second;
					temp.seq2 = itB->second;
					temp.dist = distance;
					sameSeqs.push_back(temp);
					exitedBreak = false;
					//what about precision??
					
				}else{ 
					next.seq1 = itA->second;
					next.seq2 = itB->second;
					next.dist = distance;
					exitedBreak = true;
					break;
				}
				
			}
		}

		return sameSeqs;
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "getSeqs");
		exit(1);
	}


}

//**********************************************************************************************************************

