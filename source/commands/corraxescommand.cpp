/*
 *  corraxescommand.cpp
 *  Mothur
 *
 *  Created by westcott on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "corraxescommand.h"
#include "linearalgebra.h"

//**********************************************************************************************************************
vector<string> CorrAxesCommand::setParameters(){	
	try {
		CommandParameter paxes("axes", "InputTypes", "", "", "none", "none", "none","corraxes",false,true,true); parameters.push_back(paxes);
		CommandParameter pshared("shared", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none","",false,false,true); parameters.push_back(pshared);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none","",false,false,true); parameters.push_back(prelabund);
		CommandParameter pmetadata("metadata", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none","",false,false); parameters.push_back(pmetadata);
		CommandParameter pnumaxes("numaxes", "Number", "", "3", "", "", "","",false,false); parameters.push_back(pnumaxes);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pmethod("method", "Multiple", "pearson-spearman-kendall", "pearson", "", "", "","",false,false); parameters.push_back(pmethod);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["corraxes"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CorrAxesCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The corr.axes command reads a shared, relabund or metadata file as well as an axes file and calculates the correlation coefficient.\n";
		helpString += "The corr.axes command parameters are shared, relabund, axes, metadata, groups, method, numaxes and label.  The shared, relabund or metadata and axes parameters are required.  If shared is given the relative abundance is calculated.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance level you would like used, if none is given the first distance is used.\n";
		helpString += "The method parameter allows you to select what method you would like to use. Options are pearson, spearman and kendall. Default=pearson.\n";
		helpString += "The numaxes parameter allows you to select the number of axes you would like to use. Default=3.\n";
		helpString += "The corr.axes command should be in the following format: corr.axes(axes=yourPcoaFile, shared=yourSharedFile, method=yourMethod).\n";
		helpString += "Example corr.axes(axes=genus.pool.thetayc.genus.lt.pcoa, shared=genus.pool.shared, method=kendall).\n";
		helpString += "The corr.axes command outputs a .corr.axes file.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string CorrAxesCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "corraxes") {  pattern = "[filename],[tag],corr.axes"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "CorrAxesCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
CorrAxesCommand::CorrAxesCommand(string option) : Command()  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			axesfile = validParameter.validFile(parameters, "axes");
			if (axesfile == "not open") { abort = true; }
			else if (axesfile == "not found") { axesfile = ""; m->mothurOut("axes is a required parameter for the corr.axes command.\n");  abort = true;  }	
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; current->setRelAbundFile(relabundfile); }
			
			metadatafile = validParameter.validFile(parameters, "metadata");
			if (metadatafile == "not open") { abort = true; }
			else if (metadatafile == "not found") { metadatafile = ""; }
			else { inputFileName = metadatafile; }
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "";  pickedGroups = false;  }
			else { 
				pickedGroups = true;
				util.splitAtDash(groups, Groups);
                    if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }	
			}			
			
					if (outputdir == ""){    outputdir = util.hasPath(inputFileName);	}
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile.\n");  label=""; }	
			
			if ((relabundfile == "") && (sharedfile == "") && (metadatafile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then relabund
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { inputFileName = sharedfile; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 
					relabundfile = current->getRelAbundFile(); 
					if (relabundfile != "") { inputFileName = relabundfile;  m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter.\n");  }
					else { 
						m->mothurOut("You must provide either a shared, relabund, or metadata file.\n");  abort = true; 
					}
				}
			}	
			
			if (metadatafile != "") {
				if ((relabundfile != "") || (sharedfile != "")) { m->mothurOut("You may only use one of the following : shared, relabund or metadata file.\n");  abort = true;  }
			}else {
				if ((relabundfile != "") && (sharedfile != "")) { m->mothurOut("You may only use one of the following : shared, relabund or metadata file.\n");  abort = true;  }
			}
			string temp;
			temp = validParameter.valid(parameters, "numaxes");		if (temp == "not found"){	temp = "3";				}
			util.mothurConvert(temp, numaxes); 
			
			method = validParameter.valid(parameters, "method");		if (method == "not found"){	method = "pearson";		}
			
			if ((method != "pearson") && (method != "spearman") && (method != "kendall")) { m->mothurOut(method + " is not a valid method. Valid methods are pearson, spearman, and kendall.\n");  abort = true; }
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "CorrAxesCommand");		
		exit(1);
	}
}
//**********************************************************************************************************************

int CorrAxesCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		/*************************************************************************************/
		// use smart distancing to get right sharedRabund and convert to relabund if needed  //
		/************************************************************************************/
		if (sharedfile != "") {  
			InputData* input = new InputData(sharedfile, "sharedfile", Groups);
			getSharedFloat(input); 
			delete input;
			
            if (m->getControl_pressed()) {  delete lookupFloat; return 0; }
			if (lookupFloat == nullptr) { m->mothurOut("[ERROR] reading relabund file.\n");  return 0; }
			
		}else if (relabundfile != "") { 
			InputData* input = new InputData(relabundfile, "relabund", Groups);
			getSharedFloat(input); 
			delete input;
			
			if (m->getControl_pressed()) {  delete lookupFloat; return 0; }
			if (lookupFloat == nullptr) { m->mothurOut("[ERROR] reading relabund file.\n");  return 0; }
			
		}else if (metadatafile != "") { 
			getMetadata();  //reads metadata file and store in lookupFloat, saves column headings in metadataLabels for later
			if (m->getControl_pressed()) {  delete lookupFloat; return 0; }
			if (lookupFloat == nullptr) { m->mothurOut("[ERROR] reading metadata file.\n");  return 0; }
		}else {	m->mothurOut("[ERROR]: no file given.\n");  return 0; }
		
		if (m->getControl_pressed()) {  delete lookupFloat; return 0; }
		
		//this is for a sanity check to make sure the axes file and shared file match
        vector<string> lookupGroups = lookupFloat->getNamesGroups();
        for (int i = 0; i < lookupGroups.size(); i++) { names.insert(lookupGroups[i]); }
		
		/*************************************************************************************/
		// read axes file and check for file mismatches with shared or relabund file         //
		/************************************************************************************/
		
		//read axes file
		map<string, vector<float> > axes = readAxes();
		
		if (m->getControl_pressed()) {  delete lookupFloat; return 0; }
		
		//sanity check, the read only adds groups that are in the shared or relabund file, but we want to make sure the axes file isn't missing anyone
		if (axes.size() != lookupGroups.size()) {
			map<string, vector<float> >::iterator it;
			for (int i = 0; i < lookupGroups.size(); i++) {
				it = axes.find(lookupGroups[i]);
				if (it == axes.end()) { m->mothurOut(lookupGroups[i] + " is in your shared of relabund file but not in your axes file, please correct.\n");  }
			}
			m->setControl_pressed(true);
		}
		
		if (m->getControl_pressed()) {  delete lookupFloat; return 0; }
		
		/*************************************************************************************/
		// calc the r values																//
		/************************************************************************************/
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[tag]"] = method;
		string outputFileName = getOutputFileName("corraxes", variables);
		outputNames.push_back(outputFileName); outputTypes["corraxes"].push_back(outputFileName);	
		ofstream out;
		util.openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//output headings
		if (metadatafile == "") {  out << "OTU";	}
		else {  out << "Feature";						}

		for (int i = 0; i < numaxes; i++) { out << '\t' << "axis" << (i+1) << "\tp-value"; }
		out << "\tlength" << endl;
		
		if (method == "pearson")		{  calcPearson(axes, out);	}
		else if (method == "spearman")	{  calcSpearman(axes, out); }
		else if (method == "kendall")	{  calcKendall(axes, out);	}
		else { m->mothurOut("[ERROR]: Invalid method.\n");  }
		
		out.close();
		delete lookupFloat;
		
		if (m->getControl_pressed()) {  return 0; }

		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "execute");	
		exit(1);
	}
}
//**********************************************************************************************************************
int CorrAxesCommand::calcPearson(map<string, vector<float> >& axes, ofstream& out) {
   try {
	   
       LinearAlgebra linear;
       vector<string> lookupGroups = lookupFloat->getNamesGroups();
       
	   //find average of each axis - X
	   vector<float> averageAxes; averageAxes.resize(numaxes, 0.0);
	   for (map<string, vector<float> >::iterator it = axes.begin(); it != axes.end(); it++) {
		   vector<float> temp = it->second;
		   for (int i = 0; i < temp.size(); i++) {
			   averageAxes[i] += temp[i];  
		   }
	   }
	   
	   for (int i = 0; i < averageAxes.size(); i++) {  averageAxes[i] = averageAxes[i] / (float) axes.size(); }
	   
       vector<string> currentLabels = lookupFloat->getOTUNames();
	   //for each otu
	   for (int i = 0; i < lookupFloat->getNumBins(); i++) {
		   
		   if (metadatafile == "") {  out << currentLabels[i];	}
		   else {  out << metadataLabels[i];		}
		   		   
		   //find the averages this otu - Y
		   float sumOtu = lookupFloat->getOTUTotal(i);
		   float Ybar = sumOtu / (float) lookupFloat->size();
		   
		   vector<float> rValues(averageAxes.size());

		   //find r value for each axis
		   for (int k = 0; k < averageAxes.size(); k++) {
			   
			   double r = 0.0;
			   double numerator = 0.0;
			   double denomTerm1 = 0.0;
			   double denomTerm2 = 0.0;
			   for (int j = 0; j < lookupFloat->size(); j++) {
				   float Yi = lookupFloat->get(i, lookupGroups[j]);
				   float Xi = axes[lookupGroups[j]][k];
				   
				   numerator += ((Xi - averageAxes[k]) * (Yi - Ybar));
				   denomTerm1 += ((Xi - averageAxes[k]) * (Xi - averageAxes[k]));
				   denomTerm2 += ((Yi - Ybar) * (Yi - Ybar));
			   }
			   
			   double denom = (sqrt(denomTerm1) * sqrt(denomTerm2));
			   
			   r = numerator / denom;
               
               if (isnan(r) || isinf(r)) { r = 0.0; }
               
			   rValues[k] = r;
			   out << '\t' << r; 
               
               double sig = linear.calcPearsonSig(lookupFloat->size(), r);
               
               out << '\t' << sig;
		   }
		   
		   double sum = 0;
		   for(int k=0;k<rValues.size();k++){
			   sum += rValues[k] * rValues[k];
		   }
		   out << '\t' << sqrt(sum) << endl;
	   }
	   	   
	   return 0;
   }
   catch(exception& e) {
	   m->errorOut(e, "CorrAxesCommand", "calcPearson");
	   exit(1);
   }
}
//**********************************************************************************************************************
int CorrAxesCommand::calcSpearman(map<string, vector<float> >& axes, ofstream& out) {
	try {
		
        LinearAlgebra linear;
        vector<string> lookupGroups = lookupFloat->getNamesGroups();
        vector<double> sf; 
        
		//format data
		vector< map<float, int> > tableX; tableX.resize(numaxes);
		map<float, int>::iterator itTable;
		vector< vector<spearmanRank> > scores; scores.resize(numaxes);
		for (map<string, vector<float> >::iterator it = axes.begin(); it != axes.end(); it++) {
			vector<float> temp = it->second;
			for (int i = 0; i < temp.size(); i++) {
				spearmanRank member(it->first, temp[i]);
				scores[i].push_back(member);  
				
				//count number of repeats
				itTable = tableX[i].find(temp[i]);
				if (itTable == tableX[i].end()) { 
					tableX[i][temp[i]] = 1;
				}else {
					tableX[i][temp[i]]++;
				}
			}
		}
		
		//calc LX
		//for each axis
		vector<double> Lx; Lx.resize(numaxes, 0.0);
		for (int i = 0; i < numaxes; i++) {
			for (itTable = tableX[i].begin(); itTable != tableX[i].end(); itTable++) {
				double tx = (double) itTable->second;
				Lx[i] += ((pow(tx, 3.0) - tx) / 12.0);
			}
		}
		
		//sort each axis
		for (int i = 0; i < numaxes; i++) {  sort(scores[i].begin(), scores[i].end(), compareSpearman); }
		
		//find ranks of xi in each axis
		map<string, vector<float> > rankAxes;
		for (int i = 0; i < numaxes; i++) {
			
			vector<spearmanRank> ties;
			int rankTotal = 0;
            double sfTemp = 0.0;
			for (int j = 0; j < scores[i].size(); j++) {
				rankTotal += (j+1);
				ties.push_back(scores[i][j]);
				
				if (j != (scores[i].size()-1)) { // you are not the last so you can look ahead
					if (scores[i][j].score != scores[i][j+1].score) { // you are done with ties, rank them and continue

						for (int k = 0; k < ties.size(); k++) {
							float thisrank = rankTotal / (float) ties.size();
  							rankAxes[ties[k].name].push_back(thisrank);
						}
                        int t = ties.size();
                        sfTemp += (t*t*t-t);
						ties.clear();
						rankTotal = 0;
					}
				}else { // you are the last one
					
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						rankAxes[ties[k].name].push_back(thisrank);
						
					}
				}
			}
            sf.push_back(sfTemp);
		}
		
		vector<string> currentLabels = lookupFloat->getOTUNames();
        
		//for each otu
		for (int i = 0; i < lookupFloat->getNumBins(); i++) {
			
			if (metadatafile == "") {  out << currentLabels[i];	}
			else {  out << metadataLabels[i];		}
			
			//find the ranks of this otu - Y
			vector<spearmanRank> otuScores;
			map<float, int> tableY;
			for (int j = 0; j < lookupFloat->size(); j++) {
				spearmanRank member(lookupGroups[j], lookupFloat->get(i, lookupGroups[j]));
				otuScores.push_back(member);
				
				itTable = tableY.find(member.score);
				if (itTable == tableY.end()) { 
					tableY[member.score] = 1;
				}else {
					tableY[member.score]++;
				}
				
			}
			
			//calc Ly
			double Ly = 0.0;
			for (itTable = tableY.begin(); itTable != tableY.end(); itTable++) {
				double ty = (double) itTable->second;
				Ly += ((pow(ty, 3.0) - ty) / 12.0);
			}
			
			sort(otuScores.begin(), otuScores.end(), compareSpearman);
			
            double sg = 0.0;
			map<string, float> rankOtus;
			vector<spearmanRank> ties;
			int rankTotal = 0;
			for (int j = 0; j < otuScores.size(); j++) {
				rankTotal += (j+1);
				ties.push_back(otuScores[j]);
				
				if (j != (otuScores.size()-1)) { // you are not the last so you can look ahead
					if (otuScores[j].score != otuScores[j+1].score) { // you are done with ties, rank them and continue
						
						for (int k = 0; k < ties.size(); k++) {
							float thisrank = rankTotal / (float) ties.size();
  							rankOtus[ties[k].name] = thisrank;
						}
                        int t = ties.size();
                        sg += (t*t*t-t);
						ties.clear();
						rankTotal = 0;
					}
				}else { // you are the last one
					
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						rankOtus[ties[k].name] = thisrank;
					}
				}
			}
			vector<double> pValues(numaxes);	

			//calc spearman ranks for each axis for this otu
			for (int j = 0; j < numaxes; j++) {
				
				double di = 0.0;
				for (int k = 0; k < lookupFloat->size(); k++) {
					
					float xi = rankAxes[lookupGroups[k]][j];
					float yi = rankOtus[lookupGroups[k]];
					
					di += ((xi - yi) * (xi - yi));
				}
				
				double p = 0.0;
				
				double n = (double) lookupFloat->size();
				
				double SX2 = ((pow(n, 3.0) - n) / 12.0) - Lx[j];
				double SY2 = ((pow(n, 3.0) - n) / 12.0) - Ly;
				
				p = (SX2 + SY2 - di) / (2.0 * sqrt((SX2*SY2)));
				
                if (isnan(p) || isinf(p)) { p = 0.0; }
                
				out  << '\t' << p;
				
				pValues[j] = p;
                
                double sig = linear.calcSpearmanSig(n, sf[j], sg, di);            
                out  << '\t' << sig;
                
			}

			double sum = 0;
			for(int k=0;k<numaxes;k++){
				sum += pValues[k] * pValues[k];
			}
			out << '\t' << sqrt(sum) << endl;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "calcSpearman");
		exit(1);
	}
}
//**********************************************************************************************************************
int CorrAxesCommand::calcKendall(map<string, vector<float> >& axes, ofstream& out) {
	try {
		
        LinearAlgebra linear;
        vector<string> lookupGroups = lookupFloat->getNamesGroups();
        
		//format data
		vector< vector<spearmanRank> > scores; scores.resize(numaxes);
		for (map<string, vector<float> >::iterator it = axes.begin(); it != axes.end(); it++) {
			vector<float> temp = it->second;
			for (int i = 0; i < temp.size(); i++) {
				spearmanRank member(it->first, temp[i]);
				scores[i].push_back(member);  
			}
		}
		
		//sort each axis
		for (int i = 0; i < numaxes; i++) {  sort(scores[i].begin(), scores[i].end(), compareSpearman); }
		
		//convert scores to ranks of xi in each axis
		for (int i = 0; i < numaxes; i++) {
			
			vector<spearmanRank*> ties;
			int rankTotal = 0;
			for (int j = 0; j < scores[i].size(); j++) {
				rankTotal += (j+1);
				ties.push_back(&(scores[i][j]));
				
				if (j != scores[i].size()-1) { // you are not the last so you can look ahead
					if (scores[i][j].score != scores[i][j+1].score) { // you are done with ties, rank them and continue
						for (int k = 0; k < ties.size(); k++) {
							float thisrank = rankTotal / (float) ties.size();
  							(*ties[k]).score = thisrank;
						}
						ties.clear();
						rankTotal = 0;
					}
				}else { // you are the last one
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						(*ties[k]).score = thisrank;
					}
				}
			}
		}
		
        vector<string> currentLabels = lookupFloat->getOTUNames();
        
		//for each otu
		for (int i = 0; i < lookupFloat->getNumBins(); i++) {
		
			if (metadatafile == "") {  out << currentLabels[i];	}
			else {  out << metadataLabels[i];		}
			
			//find the ranks of this otu - Y
			vector<spearmanRank> otuScores;
			for (int j = 0; j < lookupFloat->size(); j++) {
				spearmanRank member(lookupGroups[j], lookupFloat->get(i, lookupGroups[j]));
				otuScores.push_back(member);
			}
						
			sort(otuScores.begin(), otuScores.end(), compareSpearman);
			
			map<string, float> rankOtus;
			vector<spearmanRank> ties;
			int rankTotal = 0;
			for (int j = 0; j < otuScores.size(); j++) {
				rankTotal += (j+1);
				ties.push_back(otuScores[j]);
				
				if (j != otuScores.size()-1) { // you are not the last so you can look ahead
					if (otuScores[j].score != otuScores[j+1].score) { // you are done with ties, rank them and continue
						for (int k = 0; k < ties.size(); k++) {
							float thisrank = rankTotal / (float) ties.size();
  							rankOtus[ties[k].name] = thisrank;
						}
						ties.clear();
						rankTotal = 0;
					}
				}else { // you are the last one
					for (int k = 0; k < ties.size(); k++) {
						float thisrank = rankTotal / (float) ties.size();
						rankOtus[ties[k].name] = thisrank;
					}
				}
			}
			
			
			vector<double> pValues(numaxes);
			
			//calc spearman ranks for each axis for this otu
			for (int j = 0; j < numaxes; j++) {
			
				int numCoor = 0;
				int numDisCoor = 0;
				
				vector<spearmanRank> otus; 
				vector<spearmanRank> otusTemp;
				for (int l = 0; l < scores[j].size(); l++) {   
					spearmanRank member(scores[j][l].name, rankOtus[scores[j][l].name]);
					otus.push_back(member);
				}
				
				int count = 0;
				for (int l = 0; l < scores[j].size(); l++) {
					
					int numWithHigherRank = 0;
					int numWithLowerRank = 0;
					float thisrank = otus[l].score;
					
					for (int u = l+1; u < scores[j].size(); u++) {
						if (otus[u].score > thisrank) { numWithHigherRank++; }
						else if (otus[u].score < thisrank) { numWithLowerRank++; }
						count++;
					}
					
					numCoor += numWithHigherRank;
					numDisCoor += numWithLowerRank;
				}
				
				double p = (numCoor - numDisCoor) / (float) count;
                 if (isnan(p) || isinf(p)) { p = 0.0; }
                
				out << '\t' << p;
				pValues[j] = p;
                
                double sig = linear.calcKendallSig(scores[j].size(), p);
                
                out << '\t' << sig;
			}
			
			double sum = 0;
			for(int k=0;k<numaxes;k++){
				sum += pValues[k] * pValues[k];
			}
			out << '\t' << sqrt(sum) << endl;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "calcKendall");
		exit(1);
	}
}
//**********************************************************************************************************************
int CorrAxesCommand::getSharedFloat(InputData* input){
	try {
		lookupFloat = input->getSharedRAbundFloatVectors();
		string lastLabel = lookupFloat->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookupFloat != nullptr) && (userLabels.size() != 0)) {
			
			if (m->getControl_pressed()) {  return 0;  }
			
			if(labels.count(lookupFloat->getLabel()) == 1){
				processedLabels.insert(lookupFloat->getLabel());
				userLabels.erase(lookupFloat->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(lookupFloat->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookupFloat->getLabel();
				
				delete lookupFloat;
				lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
				
				processedLabels.insert(lookupFloat->getLabel());
				userLabels.erase(lookupFloat->getLabel());
				
				//restore real lastlabel to save below
				lookupFloat->setLabels(saveLabel);
				break;
			}
			
			lastLabel = lookupFloat->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete lookupFloat;
			lookupFloat = input->getSharedRAbundFloatVectors();
		}
		
		
		if (m->getControl_pressed()) { return 0;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
            delete lookupFloat;
			lookupFloat = input->getSharedRAbundFloatVectors();
		}	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "getSharedFloat");	
		exit(1);
	}
}
/*****************************************************************/
map<string, vector<float> > CorrAxesCommand::readAxes(){
	try {
		map<string, vector<float> > axes;
		
		ifstream in; util.openInputFile(axesfile, in);
		
		string headerLine = util.getline(in); gobble(in);
		
		//count the number of axis you are reading
		bool done = false;
		int count = 0;
		while (!done) {
			int pos = headerLine.find("axis");
			if (pos != string::npos) {
				count++;
				headerLine = headerLine.substr(pos+4);
			}else { done = true; }
		}
		
		if (numaxes > count) { m->mothurOut("You requested " + toString(numaxes) + " axes, but your file only includes " + toString(count) + ". Using " + toString(count) + ".\n");  numaxes = count; }
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return axes; }
			
			string group = "";
			in >> group; gobble(in);
			
			vector<float> thisGroupsAxes;
			for (int i = 0; i < count; i++) {
				float temp = 0.0;
				in >> temp; 
				
				//only save the axis we want
				if (i < numaxes) {  thisGroupsAxes.push_back(temp); }
			}
			
			//save group if its one the user selected
			if (names.count(group) != 0) {
				map<string, vector<float> >::iterator it = axes.find(group);
				
				if (it == axes.end()) {
					axes[group] = thisGroupsAxes;
                    
				}else { m->mothurOut(group + " is already in your axes file, using first definition.\n");  }
			}
			
			gobble(in);
		}
		in.close();
		
		return axes;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "readAxes");	
		exit(1);
	}
}
/*****************************************************************/
int CorrAxesCommand::getMetadata(){
	try {
		vector<string> groupNames;
		
		ifstream in; util.openInputFile(metadatafile, in);
		
		string headerLine = util.getline(in); gobble(in);
		vector<string> pieces = util.splitWhiteSpace(headerLine);
		
		//save names of columns you are reading
		for (int i = 1; i < pieces.size(); i++) { metadataLabels.push_back(pieces[i]); }
			
		//read rest of file
        lookupFloat = new SharedRAbundFloatVectors();
        lookupFloat->setOTUNames(metadataLabels);
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return 0; }
			
            string metadataLine = util.getline(in); gobble(in);
            vector<string> metaPieces = util.splitWhiteSpace(metadataLine);
			string group = metaPieces[0];
			
			groupNames.push_back(group);
            
            if (pieces.size() != metaPieces.size()) { //sanity check
                m->mothurOut("[ERROR]: Your metadata file seems to be mismatched. I read " + toString(pieces.size()) + " column headers, but found " + toString(metaPieces.size()) + " columns for the " + group + " group. Quitting.\n"); m->setControl_pressed(true); in.close(); return 0;
            }
				
			SharedRAbundFloatVector* tempLookup = new SharedRAbundFloatVector();
			tempLookup->setLabel("1");
            tempLookup->setGroup(group);
        
			for (int i = 1; i < metaPieces.size(); i++) {
				float temp = 0.0;
                util.mothurConvert(metaPieces[i], temp);
				tempLookup->push_back(temp);
			}
            lookupFloat->push_back(tempLookup);
		}
		in.close();
		
        lookupFloat->setLabels("1");
        lookupFloat->eliminateZeroOTUS();
        metadataLabels = lookupFloat->getOTUNames();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CorrAxesCommand", "getMetadata");	
		exit(1);
	}
}
/*****************************************************************/






