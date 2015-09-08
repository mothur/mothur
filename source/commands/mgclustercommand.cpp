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
vector<string> MGClusterCommand::setParameters(){	
	try {
		CommandParameter pblast("blast", "InputTypes", "", "", "none", "none", "none","list",false,true,true); parameters.push_back(pblast);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName","rabund-sabund",false,false,true); parameters.push_back(pname);
		CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter plength("length", "Number", "", "5", "", "", "","",false,false); parameters.push_back(plength);
		CommandParameter ppenalty("penalty", "Number", "", "0.10", "", "", "","",false,false); parameters.push_back(ppenalty);
		CommandParameter pcutoff("cutoff", "Number", "", "0.70", "", "", "","",false,false,true); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average", "average", "", "", "","",false,false); parameters.push_back(pmethod);
		CommandParameter phard("hard", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(phard);
		CommandParameter pmin("min", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pmin);
		CommandParameter pmerge("merge", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pmerge);
        CommandParameter padjust("adjust", "String", "", "F", "", "", "","",false,false); parameters.push_back(padjust);
		CommandParameter phcluster("hcluster", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(phcluster);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MGClusterCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The mgcluster command parameter options are blast, name, cutoff, precision, hard,  method, merge, min, length, penalty, adjust and hcluster. The blast parameter is required.\n";
		helpString += "The mgcluster command reads a blast and name file and clusters the sequences into OPF units similar to the OTUs.\n";
		helpString += "This command outputs a .list, .rabund and .sabund file that can be used with mothur other commands to estimate richness.\n";
		helpString += "The cutoff parameter is used to specify the maximum distance you would like to cluster to. The default is 0.70.\n";
		helpString += "The precision parameter's default value is 100. \n";
		helpString += "The acceptable mgcluster methods are furthest, nearest and average.  If no method is provided then average is assumed.\n";	
		helpString += "The min parameter allows you to specify is you want the minimum or maximum blast score ratio used in calculating the distance. The default is true, meaning you want the minimum.\n";
		helpString += "The length parameter is used to specify the minimum overlap required.  The default is 5.\n";
        helpString += "The adjust parameter is used to handle missing distances.  If you set a cutoff, adjust=f by default.  If not, adjust=t by default. Adjust=f, means ignore missing distances and adjust cutoff as needed with the average neighbor method.  Adjust=t, will treat missing distances as 1.0. You can also set the value the missing distances should be set to, adjust=0.5 would give missing distances a value of 0.5.\n";
		helpString += "The penalty parameter is used to adjust the error rate.  The default is 0.10.\n";
		helpString += "The merge parameter allows you to shut off merging based on overlaps and just cluster.  By default merge is true, meaning you want to merge.\n";
		helpString += "The hcluster parameter allows you to use the hcluster algorithm when clustering.  This may be necessary if your file is too large to fit into RAM. The default is false.\n";
		helpString += "The mgcluster command should be in the following format: \n";
		helpString += "mgcluster(blast=yourBlastfile, name=yourNameFile, cutoff=yourCutOff).\n";
		helpString += "Note: No spaces between parameter labels (i.e. balst), '=' and parameters (i.e.yourBlastfile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MGClusterCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "list") {  pattern = "[filename],[clustertag],list-[filename],[clustertag],[tag2],list"; } 
        else if (type == "rabund") {  pattern = "[filename],[clustertag],rabund"; } 
        else if (type == "sabund") {  pattern = "[filename],[clustertag],sabund"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MGClusterCommand", "getOutputPattern");
        exit(1);
    }
}
//*******************************************************************************************************************
MGClusterCommand::MGClusterCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "MGClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
MGClusterCommand::MGClusterCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("blast");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["blast"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			blastfile = validParameter.validFile(parameters, "blast", true);
			if (blastfile == "not open") { blastfile = ""; abort = true; }	
			else if (blastfile == "not found") { blastfile = ""; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(blastfile); //if user entered a file with a path then preserve it	
			}
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; }	
			else if (countfile == "not found") { countfile = ""; }
            else { m->setCountTableFile(countfile); }
            
            if (countfile != "" && namefile != "") { m->mothurOut("[ERROR]: Cannot have both a name file and count file. Please use one or the other."); m->mothurOutEndLine(); abort = true; }
			
			if ((blastfile == "")) { m->mothurOut("When executing a mgcluster command you must provide a blastfile."); m->mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			string temp;
            temp = validParameter.validFile(parameters, "precision", false);		if (temp == "not found") { temp = "100"; }
			precisionLength = temp.length();
			m->mothurConvert(temp, precision); 
			
            cutoffSet = false;
			temp = validParameter.validFile(parameters, "cutoff", false);
            if (temp == "not found") { temp = "0.70"; }
            else { cutoffSet = true;  }
			m->mothurConvert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "average"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average."); m->mothurOutEndLine(); abort = true; }

			temp = validParameter.validFile(parameters, "length", false);			if (temp == "not found") { temp = "5"; }
			m->mothurConvert(temp, length); 
			
			temp = validParameter.validFile(parameters, "penalty", false);			if (temp == "not found") { temp = "0.10"; }
			m->mothurConvert(temp, penalty); 
			
			temp = validParameter.validFile(parameters, "min", false);				if (temp == "not found") { temp = "true"; }
			minWanted = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "merge", false);			if (temp == "not found") { temp = "true"; }
			merge = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "hcluster", false);			if (temp == "not found") { temp = "false"; }
			hclusterWanted = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "T"; }
			hard = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "adjust", false);				if (temp == "not found") { if (cutoffSet) { temp = "F"; }else { temp="T"; } }
            if (m->isNumeric1(temp))    { m->mothurConvert(temp, adjust);   }
            else if (m->isTrue(temp))   { adjust = 1.0;                     }
            else                        { adjust = -1.0;                    }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "MGClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int MGClusterCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read names file
        map<string, int> counts;
		if (namefile != "") {
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
		}else if (countfile != "") {
            ct = new CountTable();
            ct->readTable(countfile, false, false);
            nameMap= new NameAssignment();
            vector<string> tempNames = ct->getNamesOfSeqs();
            for (int i = 0; i < tempNames.size(); i++) {  nameMap->push_back(tempNames[i]);  }
            counts = ct->getNameMap();
        }else{ nameMap= new NameAssignment(); }
		
		string fileroot = outputDir + m->getRootName(m->getSimpleName(blastfile));
		string tag = "";
		time_t start;
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000; 
        
		//read blastfile - creates sparsematrices for the distances and overlaps as well as a listvector
		//must remember to delete those objects here since readBlast does not
		read = new ReadBlast(blastfile, cutoff, penalty, length, minWanted, hclusterWanted);
		read->read(nameMap);
        
        list = new ListVector(nameMap->getListVector());
        RAbundVector* rabund = NULL;
        
        if(countfile != "") {
            rabund = new RAbundVector();
            createRabund(ct, list, rabund);
        }else {
            rabund = new RAbundVector(list->getRAbundVector());
        }
        
                
		//list = new ListVector(nameMap->getListVector());
		//rabund = new RAbundVector(list->getRAbundVector());
		
		if (m->control_pressed) { outputTypes.clear(); delete nameMap; delete read; delete list; delete rabund; return 0; }
		
		start = time(NULL);
		oldList = *list;
		map<string, int> Seq2Bin;
		map<string, int> oldSeq2Bin;
		
		if (method == "furthest")		{ tag = "fn";  }
		else if (method == "nearest")	{ tag = "nn";  }
		else							{ tag = "an";  }
		
        map<string, string> variables; 
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	sabundFile);
            m->openOutputFile(rabundFileName,	rabundFile);
        }
		m->openOutputFile(listFileName,	listFile);
        list->printHeaders(listFile);
		
		if (m->control_pressed) { 
			delete nameMap; delete read; delete list; delete rabund; 
			listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
			outputTypes.clear();
			return 0; 
		}
		
		double saveCutoff = cutoff;
		
		if (!hclusterWanted) {
			//get distmatrix and overlap
			SparseDistanceMatrix* distMatrix = read->getDistMatrix();
			overlapMatrix = read->getOverlapMatrix(); //already sorted by read 
			delete read;
		
			//create cluster
			if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, distMatrix, cutoff, method, adjust); }
			else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, distMatrix, cutoff, method, adjust); }
			else if(method == "average"){	cluster = new AverageLinkage(rabund, list, distMatrix, cutoff, method, adjust);	}
			cluster->setMapWanted(true);
			Seq2Bin = cluster->getSeqtoBin();
			oldSeq2Bin = Seq2Bin;
			
			if (m->control_pressed) { 
				delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster;
				listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
				outputTypes.clear();
				return 0; 
			}
            
            
			//cluster using cluster classes
			while (distMatrix->getSmallDist() < cutoff && distMatrix->getNNodes() > 0){
				
                if (m->debug) {  cout << "numNodes=" << distMatrix->getNNodes() << " smallDist = " << distMatrix->getSmallDist() << endl; }
                
				cluster->update(cutoff);
				
				if (m->control_pressed) { 
					delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster;
					listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
					outputTypes.clear();
					return 0; 
				}
				
				float dist = distMatrix->getSmallDist();
				float rndDist;
				if (hard) {
					rndDist = m->ceilDist(dist, precision); 
				}else{
					rndDist = m->roundDist(dist, precision); 
				}
				
				if(previousDist <= 0.0000 && dist != previousDist){
					oldList.setLabel("unique");
					printData(&oldList, counts);
				}
				else if(rndDist != rndPreviousDist){
					if (merge) {
						ListVector* temp = mergeOPFs(oldSeq2Bin, rndPreviousDist);
						
						if (m->control_pressed) { 
							delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster; delete temp;
							listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
							outputTypes.clear();
							return 0; 
						}
						
						temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
						printData(temp, counts);
						delete temp;
					}else{
						oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
						printData(&oldList, counts);
					}
				}
	
				previousDist = dist;
				rndPreviousDist = rndDist;
				oldList = *list;
				Seq2Bin = cluster->getSeqtoBin();
				oldSeq2Bin = Seq2Bin;
			}
			
			if(previousDist <= 0.0000){
				oldList.setLabel("unique");
				printData(&oldList, counts);
			}
			else if(rndPreviousDist<cutoff){
				if (merge) {
					ListVector* temp = mergeOPFs(oldSeq2Bin, rndPreviousDist);
					
					if (m->control_pressed) { 
							delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster; delete temp;
							listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
							outputTypes.clear();
							return 0; 
					}
					
					temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(temp, counts);
					delete temp;
				}else{
					oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(&oldList, counts);
				}
			}
			
			//free memory
			overlapMatrix.clear();
			delete distMatrix;
			delete cluster;
			
		}else { //use hcluster to cluster
			//get distmatrix and overlap
			overlapFile = read->getOverlapFile();
			distFile = read->getDistFile(); 
			delete read;
		
			//sort the distance and overlap files
			sortHclusterFiles(distFile, overlapFile);
			
			if (m->control_pressed) { 
				delete nameMap;  delete list; delete rabund; 
				listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
				outputTypes.clear();
				return 0; 
			}
		
			//create cluster
			hcluster = new HCluster(rabund, list, method, distFile, nameMap, cutoff);
			hcluster->setMapWanted(true);
			Seq2Bin = cluster->getSeqtoBin();
			oldSeq2Bin = Seq2Bin;
			
			vector<seqDist> seqs; seqs.resize(1); // to start loop
			//ifstream inHcluster;
			//m->openInputFile(distFile, inHcluster);
			
			if (m->control_pressed) { 
				delete nameMap;  delete list; delete rabund; delete hcluster;
				listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
				outputTypes.clear();
				return 0; 
			}

			while (seqs.size() != 0){
		
				seqs = hcluster->getSeqs();
				
				//to account for cutoff change in average neighbor
				if (seqs.size() != 0) {
					if (seqs[0].dist > cutoff) { break; }
				}
				
				if (m->control_pressed) { 
					delete nameMap;  delete list; delete rabund; delete hcluster;
					listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
					m->mothurRemove(distFile);
					m->mothurRemove(overlapFile);
					outputTypes.clear();
					return 0; 
				}
				
				for (int i = 0; i < seqs.size(); i++) {  //-1 means skip me
					
					if (seqs[i].seq1 != seqs[i].seq2) {
		
						cutoff = hcluster->update(seqs[i].seq1, seqs[i].seq2, seqs[i].dist);
						
						if (m->control_pressed) { 
							delete nameMap;  delete list; delete rabund; delete hcluster;
							listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
							m->mothurRemove(distFile);
							m->mothurRemove(overlapFile);
							outputTypes.clear();
							return 0; 
						}
	
						float rndDist;
						if (hard) {
							rndDist = m->ceilDist(seqs[i].dist, precision); 
						}else{
							rndDist = m->roundDist(seqs[i].dist, precision); 
						}
												
						if((previousDist <= 0.0000) && (seqs[i].dist != previousDist)){
							oldList.setLabel("unique");
							printData(&oldList, counts);
						}
						else if((rndDist != rndPreviousDist)){
							if (merge) {
								ListVector* temp = mergeOPFs(oldSeq2Bin, rndPreviousDist);
								
								if (m->control_pressed) { 
									delete nameMap;  delete list; delete rabund; delete hcluster; delete temp;
									listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
									m->mothurRemove(distFile);
									m->mothurRemove(overlapFile);
									outputTypes.clear();
									return 0; 
								}

								temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
								printData(temp, counts);
								delete temp;
							}else{
								oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
								printData(&oldList, counts);
							}
						}
						
						previousDist = seqs[i].dist;
						rndPreviousDist = rndDist;
						oldList = *list;
						Seq2Bin = cluster->getSeqtoBin();
						oldSeq2Bin = Seq2Bin;
					}
				}
			}
			//inHcluster.close();
			
			if(previousDist <= 0.0000){
				oldList.setLabel("unique");
				printData(&oldList, counts);
			}
			else if(rndPreviousDist<cutoff){
				if (merge) {
					ListVector* temp = mergeOPFs(oldSeq2Bin, rndPreviousDist);
					
					if (m->control_pressed) { 
							delete nameMap; delete list; delete rabund; delete hcluster; delete temp;
							listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
							m->mothurRemove(distFile);
							m->mothurRemove(overlapFile);
							outputTypes.clear();
							return 0; 
					}
					
					temp->setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(temp, counts);
					delete temp;
				}else{
					oldList.setLabel(toString(rndPreviousDist,  precisionLength-1));
					printData(&oldList, counts);
				}
			}
			
			delete hcluster;
			m->mothurRemove(distFile);
			m->mothurRemove(overlapFile);
		}
		
		delete list;
		delete rabund;
		listFile.close();
        if (countfile == "") {
            sabundFile.close();
            rabundFile.close();
        }
		if (m->control_pressed) { 
			delete nameMap; 
			listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); } m->mothurRemove((fileroot+ tag + ".list"));
			outputTypes.clear();
			return 0; 
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(listFileName); m->mothurOutEndLine();	outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
		if (countfile == "") {
            m->mothurOut(rabundFileName); m->mothurOutEndLine();	outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            m->mothurOut(sabundFileName); m->mothurOutEndLine();	outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
        }
		m->mothurOutEndLine();
		
		if (saveCutoff != cutoff) { 
			if (hard)	{  saveCutoff = m->ceilDist(saveCutoff, precision);	}
			else		{	saveCutoff = m->roundDist(saveCutoff, precision);  }
			
			m->mothurOut("changed cutoff to " + toString(cutoff)); m->mothurOutEndLine(); 
		}
		
		//set list file as new current listfile
		string current = "";
		itTypes = outputTypes.find("list");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
		}
		
		//set rabund file as new current rabundfile
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
		}
		
		//set sabund file as new current sabundfile
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSabundFile(current); }
		}
		
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " seconds to cluster."); m->mothurOutEndLine();
			
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void MGClusterCommand::printData(ListVector* mergedList, map<string, int>& counts){
	try {
        if (countfile != "") {
            mergedList->print(listFile, counts);
        }else { mergedList->print(listFile); }
        
        SAbundVector sabund = mergedList->getSAbundVector();
        
        if (countfile == "") {
            mergedList->getRAbundVector().print(rabundFile);
            sabund.print(sabundFile);
        }

		sabund.print(cout);
	}
	catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "printData");
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
			m->openInputFile(overlapFile, inOverlap);  
			if (inOverlap.eof()) {  done = true;  }
		}else { if (overlapMatrix.size() == 0)  {  done = true;  } } 
		
		while (!done) {
			if (m->control_pressed) { 
				if (hclusterWanted) { 	inOverlap.close();  }		
				return newList;
			}
			
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
					inOverlap >> firstName >> secondName >> overlapDistance; m->gobble(inOverlap);
					
					//commented out because we check this in readblast already
					//map<string,int>::iterator itA = nameMap->find(firstName);
					//map<string,int>::iterator itB = nameMap->find(secondName);
					//if(itA == nameMap->end()){  cerr << "AAError: Sequence '" << firstName << "' was not found in the names file, please correct\n"; exit(1);  }
					//if(itB == nameMap->end()){  cerr << "ABError: Sequence '" << secondName << "' was not found in the names file, please correct\n"; exit(1);  }
					
					//overlapNode.seq1 = itA->second;
					//overlapNode.seq2 = itB->second;
					overlapNode.seq1 = nameMap->get(firstName);
					overlapNode.seq2 = nameMap->get(secondName);
					overlapNode.dist = overlapDistance;
				}else { inOverlap.close(); break; }
			} 
		
			if (overlapNode.dist < dist) {
				//get names of seqs that overlap
				string name1 = nameMap->get(overlapNode.seq1);
				string name2 = nameMap->get(overlapNode.seq2);
			
				//use binInfo to find out if they are already in the same bin
				//map<string, int>::iterator itBin1 = binInfo.find(name1);
				//map<string, int>::iterator itBin2 = binInfo.find(name2);
				
				//if(itBin1 == binInfo.end()){  cerr << "AAError: Sequence '" << name1 << "' does not have any bin info.\n"; exit(1);  }
				//if(itBin2 == binInfo.end()){  cerr << "ABError: Sequence '" << name2 << "' does not have any bin info.\n"; exit(1);  }

				//int binKeep = itBin1->second;
				//int binRemove = itBin2->second;
				
				int binKeep = binInfo[name1];
				int binRemove = binInfo[name2];
			
				//if not merge bins and update binInfo
				if(binKeep != binRemove) {
		
					//save names in old bin
					string names = newList->get(binRemove);
		
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
		m->errorOut(e, "MGClusterCommand", "mergeOPFs");
		exit(1);
	}
}
//**********************************************************************************************************************
void MGClusterCommand::sortHclusterFiles(string unsortedDist, string unsortedOverlap) {
	try {
		//sort distFile
		string sortedDistFile = m->sortFile(unsortedDist, outputDir);
		m->mothurRemove(unsortedDist);  //delete unsorted file
		distFile = sortedDistFile;
		
		//sort overlap file
		string sortedOverlapFile = m->sortFile(unsortedOverlap, outputDir);
		m->mothurRemove(unsortedOverlap);  //delete unsorted file
		overlapFile = sortedOverlapFile;
	}
	catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "sortHclusterFiles");
		exit(1);
	}
}

//**********************************************************************************************************************

void MGClusterCommand::createRabund(CountTable*& ct, ListVector*& list, RAbundVector*& rabund){
    try {
        //vector<string> names = ct.getNamesOfSeqs();

        //for ( int i; i < ct.getNumGroups(); i++ ) {    rav.push_back( ct.getNumSeqs(names[i]) );    }
        //return rav;
        
        for(int i = 0; i < list->getNumBins(); i++) { 
           vector<string> binNames;
           string bin = list->get(i);
           m->splitAtComma(bin, binNames);
           int total = 0;
           for (int j = 0; j < binNames.size(); j++) { 
               total += ct->getNumSeqs(binNames[j]);
           }
           rabund->push_back(total);   
       }
        
        
    }
    catch(exception& e) {
		m->errorOut(e, "MGClusterCommand", "createRabund");
		exit(1);
	}
    
}

//**********************************************************************************************************************


