/*
 *  mgclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 12/11/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mgclustercommand.h"

//**********************************************************************************************************************
MGClusterCommand::MGClusterCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"blast", "method", "name", "cutoff", "precision", "length", "min", "penalty", "hcluster","merge"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			blastfile = validParameter.validFile(parameters, "blast", true);
			if (blastfile == "not open") { abort = true; }	
			else if (blastfile == "not found") { blastfile = ""; }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			
			if ((blastfile == "")) { mothurOut("When executing a mgcluster command you must provide a blastfile."); mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			string temp;
			temp = validParameter.validFile(parameters, "precision", false);		if (temp == "not found") { temp = "100"; }
			precisionLength = temp.length();
			convert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "0.70"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "furthest"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else { mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average."); mothurOutEndLine(); abort = true; }

			temp = validParameter.validFile(parameters, "length", false);			if (temp == "not found") { temp = "5"; }
			convert(temp, length); 
			
			temp = validParameter.validFile(parameters, "penalty", false);			if (temp == "not found") { temp = "0.10"; }
			convert(temp, penalty); 
			
			temp = validParameter.validFile(parameters, "min", false);				if (temp == "not found") { temp = "true"; }
			minWanted = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "merge", false);			if (temp == "not found") { temp = "true"; }
			merge = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "hcluster", false);			if (temp == "not found") { temp = "false"; }
			hclusterWanted = isTrue(temp); 
		}

	}
	catch(exception& e) {
		errorOut(e, "MGClusterCommand", "MGClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void MGClusterCommand::help(){
	try {
		mothurOut("The mgcluster command parameter options are blast, name, cutoff, precision, method, merge, min, length, penalty and hcluster. The blast parameter is required.\n");
		mothurOut("The mgcluster command reads a blast and name file and clusters the sequences into OPF units similiar to the OTUs.\n");
		mothurOut("This command outputs a .list, .rabund and .sabund file that can be used with mothur other commands to estimate richness.\n");
		mothurOut("The cutoff parameter is used to specify the maximum distance you would like to cluster to. The default is 0.70.\n");
		mothurOut("The precision parameter's default value is 100. \n");
		mothurOut("The acceptable mgcluster methods are furthest, nearest and average.  If no method is provided then furthest is assumed.\n\n");	
		mothurOut("The min parameter allows you to specify is you want the minimum or maximum blast score ratio used in calculating the distance. The default is true, meaning you want the minimum.\n");
		mothurOut("The length parameter is used to specify the minimum overlap required.  The default is 5.\n");
		mothurOut("The penalty parameter is used to adjust the error rate.  The default is 0.10.\n");
		mothurOut("The merge parameter allows you to shut off merging based on overlaps and just cluster.  By default merge is true, meaning you want to merge.\n");
		mothurOut("The hcluster parameter allows you to use the hcluster algorithm when clustering.  This may be neccessary if your file is too large to fit into RAM. The default is false.\n");
		mothurOut("The mgcluster command should be in the following format: \n");
		mothurOut("mgcluster(blast=yourBlastfile, name=yourNameFile, cutoff=yourCutOff).\n");
		mothurOut("Note: No spaces between parameter labels (i.e. balst), '=' and parameters (i.e.yourBlastfile).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "MGClusterCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
MGClusterCommand::~MGClusterCommand(){}
//**********************************************************************************************************************
int MGClusterCommand::execute(){
	try {
		
		if (abort == true) {	return 0;	}
		
		//read names file
		if (namefile != "") {
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
		}else{ nameMap= new NameAssignment(); }
		
		string fileroot = getRootName(blastfile);
		string tag = "";
		time_t start;
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		
		//read blastfile - creates sparsematrices for the distances and overlaps as well as a listvector
		//must remember to delete those objects here since readBlast does not
		read = new ReadBlast(blastfile, cutoff, penalty, length, minWanted, hclusterWanted);
		read->read(nameMap);
		
		list = new ListVector(nameMap->getListVector());
		RAbundVector* rabund = new RAbundVector(list->getRAbundVector());
		
		start = time(NULL);
		oldList = *list;
		
		if (method == "furthest")		{ tag = "fn";  }
		else if (method == "nearest")	{ tag = "nn";  }
		else							{ tag = "an";  }
		
		//open output files
		openOutputFile(fileroot+ tag + ".list",  listFile);
		openOutputFile(fileroot+ tag + ".rabund",  rabundFile);
		openOutputFile(fileroot+ tag + ".sabund",  sabundFile);
		
		if (!hclusterWanted) {
			//get distmatrix and overlap
			SparseMatrix* distMatrix = read->getDistMatrix();
			overlapMatrix = read->getOverlapMatrix(); //already sorted by read 
			delete read;
		
			//create cluster
			if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, distMatrix, cutoff); }
			else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, distMatrix, cutoff); }
			else if(method == "average"){	cluster = new AverageLinkage(rabund, list, distMatrix, cutoff);	}
			cluster->setMapWanted(true);
			
			//cluster using cluster classes
			while (distMatrix->getSmallDist() < cutoff && distMatrix->getNNodes() > 0){
				
				cluster->update();
				float dist = distMatrix->getSmallDist();
				float rndDist = roundDist(dist, precision);
				
				if(previousDist <= 0.0000 && dist != previousDist){
					oldList.setLabel("unique");
					printData(&oldList);
				}
				else if(rndDist != rndPreviousDist){
					if (merge) {
						map<string, int> seq2Bin = cluster->getSeqtoBin();
						ListVector* temp = mergeOPFs(seq2Bin, rndPreviousDist);
						temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
						printData(temp);
						delete temp;
					}else{
						oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
						printData(&oldList);
					}
				}
				
				previousDist = dist;
				rndPreviousDist = rndDist;
				oldList = *list;
			}
			
			if(previousDist <= 0.0000){
				oldList.setLabel("unique");
				printData(&oldList);
			}
			else if(rndPreviousDist<cutoff){
				if (merge) {
					map<string, int> seq2Bin = cluster->getSeqtoBin();
					ListVector* temp = mergeOPFs(seq2Bin, rndPreviousDist);
					temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(temp);
					delete temp;
				}else{
					oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(&oldList);
				}
			}
			
			//free memory
			overlapMatrix.clear();
			delete distMatrix;
			delete cluster;
			
		}else { //use hcluster to cluster
			tag = "fn";
			//get distmatrix and overlap
			overlapFile = read->getOverlapFile();
			distFile = read->getDistFile(); 
			delete read;
		
			//sort the distance and overlap files
			sortHclusterFiles(distFile, overlapFile);
		
			//create cluster
			hcluster = new HCluster(rabund, list, method, distFile, nameMap, cutoff);
			hcluster->setMapWanted(true);
			
			vector<seqDist> seqs; seqs.resize(1); // to start loop
			//ifstream inHcluster;
			//openInputFile(distFile, inHcluster);

			while (seqs.size() != 0){
		
				seqs = hcluster->getSeqs();
				
				for (int i = 0; i < seqs.size(); i++) {  //-1 means skip me
					
					if (seqs[i].seq1 != seqs[i].seq2) {
		
						hcluster->update(seqs[i].seq1, seqs[i].seq2, seqs[i].dist);
	
						float rndDist = roundDist(seqs[i].dist, precision);
												
						if((previousDist <= 0.0000) && (seqs[i].dist != previousDist)){
							oldList.setLabel("unique");
							printData(&oldList);
						}
						else if((rndDist != rndPreviousDist)){
							if (merge) {
								map<string, int> seq2Bin = hcluster->getSeqtoBin();
								ListVector* temp = mergeOPFs(seq2Bin, rndPreviousDist);
								temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
								printData(temp);
								delete temp;
							}else{
								oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
								printData(&oldList);
							}
						}
						
						previousDist = seqs[i].dist;
						rndPreviousDist = rndDist;
						oldList = *list;
					}
				}
			}
			//inHcluster.close();
			
			if(previousDist <= 0.0000){
				oldList.setLabel("unique");
				printData(&oldList);
			}
			else if(rndPreviousDist<cutoff){
				if (merge) {
					map<string, int> seq2Bin = hcluster->getSeqtoBin();
					ListVector* temp = mergeOPFs(seq2Bin, rndPreviousDist);
					temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(temp);
					delete temp;
				}else{
					oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(&oldList);
				}
			}
			
			delete hcluster;
			remove(distFile.c_str());
			remove(overlapFile.c_str());
		}
		
		delete list; 
		delete rabund;
		listFile.close();
		sabundFile.close();
		rabundFile.close();
	
		globaldata->setListFile(fileroot+ tag + ".list");
		globaldata->setFormat("list");
			
		mothurOut("It took " + toString(time(NULL) - start) + " seconds to cluster."); mothurOutEndLine();
			
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MGClusterCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void MGClusterCommand::printData(ListVector* mergedList){
	try {
		mergedList->print(listFile);
		mergedList->getRAbundVector().print(rabundFile);
		
		SAbundVector sabund = mergedList->getSAbundVector();

		sabund.print(cout);
		sabund.print(sabundFile);
	}
	catch(exception& e) {
		errorOut(e, "MGClusterCommand", "printData");
		exit(1);
	}
}
//**********************************************************************************************************************
//this merging is just at the reporting level, after this info is printed to the file it is gone and does not effect the datastructures
//that are used to cluster by distance.  this is done so that the overlapping data does not have more influenece than the distance data.
ListVector* MGClusterCommand::mergeOPFs(map<string, int> binInfo, float dist){
	try {
		//create new listvector so you don't overwrite the clustering
		ListVector* newList = new ListVector(oldList);
		bool done = false;
		ifstream inOverlap;
		int count = 0;
		
		if (hclusterWanted) {  
			openInputFile(overlapFile, inOverlap);  
			if (inOverlap.eof()) {  done = true;  }
		}else { if (overlapMatrix.size() == 0)  {  done = true;  } } 
		
		while (!done) {
			
			//get next overlap
			seqDist overlapNode;
			if (!hclusterWanted) {  
				if (count < overlapMatrix.size()) { //do we have another node in the matrix
					overlapNode = overlapMatrix[count];
					count++;
				}else { break; }
			}else { 
				if (!inOverlap.eof()) {
					string firstName, secondName;
					float overlapDistance;
					inOverlap >> firstName >> secondName >> overlapDistance; gobble(inOverlap);
					
					map<string,int>::iterator itA = nameMap->find(firstName);
					map<string,int>::iterator itB = nameMap->find(secondName);
					if(itA == nameMap->end()){  cerr << "AAError: Sequence '" << firstName << "' was not found in the names file, please correct\n"; exit(1);  }
					if(itB == nameMap->end()){  cerr << "ABError: Sequence '" << secondName << "' was not found in the names file, please correct\n"; exit(1);  }
					
					overlapNode.seq1 = itA->second;
					overlapNode.seq2 = itB->second;
					overlapNode.dist = overlapDistance;
				}else { inOverlap.close(); break; }
			} 
		
			if (overlapNode.dist < dist) {
				//get names of seqs that overlap
				string name1 = nameMap->get(overlapNode.seq1);
				string name2 = nameMap->get(overlapNode.seq2);
				
				//use binInfo to find out if they are already in the same bin
				int binKeep = binInfo[name1];
				int binRemove = binInfo[name2];
				
				//if not merge bins and update binInfo
				if(binKeep != binRemove) {
					//save names in old bin
					string names = list->get(binRemove);
					
					//merge bins into name1s bin
					newList->set(binKeep, newList->get(binRemove)+','+newList->get(binKeep));
					newList->set(binRemove, "");	
					
					//update binInfo
					while (names.find_first_of(',') != -1) { 
						//get name from bin
						string name = names.substr(0,names.find_first_of(','));
						//save name and bin number
						binInfo[name] = binKeep;
						names = names.substr(names.find_first_of(',')+1, names.length());
					}
					
					//get last name
					binInfo[names] = binKeep;
				}
				
			}else { done = true; }
		}
		
		//return listvector
		return newList;
				
	}
	catch(exception& e) {
		errorOut(e, "MGClusterCommand", "mergeOPFs");
		exit(1);
	}
}
//**********************************************************************************************************************
void MGClusterCommand::sortHclusterFiles(string unsortedDist, string unsortedOverlap) {
	try {
		//sort distFile
		string sortedDistFile = sortFile(unsortedDist);
		remove(unsortedDist.c_str());  //delete unsorted file
		distFile = sortedDistFile;
		
		//sort overlap file
		string sortedOverlapFile = sortFile(unsortedOverlap);
		remove(unsortedOverlap.c_str());  //delete unsorted file
		overlapFile = sortedOverlapFile;
	}
	catch(exception& e) {
		errorOut(e, "MGClusterCommand", "sortHclusterFiles");
		exit(1);
	}
}

//**********************************************************************************************************************





