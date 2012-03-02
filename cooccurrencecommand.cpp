/*
 *  cooccurrencecommand.cpp
 *  Mothur
 *
 *  Created by kiverson on 1/2/12.
 *  Copyright 2012 Schloss Lab. All rights reserved.
 *
 */

#include "cooccurrencecommand.h"

//**********************************************************************************************************************
vector<string> CooccurrenceCommand::setParameters() {	
	try { 
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);		
		CommandParameter pmetric("metric", "Multiple", "cscore-checker-combo-vratio", "cscore", "", "", "",false,false); parameters.push_back(pmetric);
		CommandParameter pmatrix("matrixmodel", "Multiple", "sim1-sim2-sim3-sim4-sim5-sim6-sim7-sim8-sim9", "sim2", "", "", "",false,false); parameters.push_back(pmatrix);
        CommandParameter pruns("iters", "Number", "", "1000", "", "", "",false,false); parameters.push_back(pruns);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
        CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CooccurrenceCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CooccurrenceCommand::getHelpString(){	
	try {
		string helpString = "help!";

		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CooccurrenceCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
CooccurrenceCommand::CooccurrenceCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
        vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;

	}
	catch(exception& e) {
		m->errorOut(e, "CooccurrenceCommand", "CooccurrenceCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
CooccurrenceCommand::CooccurrenceCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
				
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}

			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
		
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;
		
	        //check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//get shared file
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(sharedfile);		}

			
			metric = validParameter.validFile(parameters, "metric", false);				if (metric == "not found") { metric = "cscore"; }
			
			if ((metric != "cscore") && (metric != "checker") && (metric != "combo") && (metric != "vratio")) {
				m->mothurOut("[ERROR]: " + metric + " is not a valid metric option for the cooccurrence command. Choices are cscore, checker, combo, vratio."); m->mothurOutEndLine(); abort = true; 
			}
			
			matrix = validParameter.validFile(parameters, "matrix", false);				if (matrix == "not found") { matrix = "sim2"; }
			
			if ((matrix != "sim1") && (matrix != "sim2") && (matrix != "sim3") && (matrix != "sim4") && (matrix != "sim5" ) && (matrix != "sim6" ) && (matrix != "sim7" ) && (matrix != "sim8" ) && (matrix != "sim9" )) {
				m->mothurOut("[ERROR]: " + matrix + " is not a valid matrix option for the cooccurrence command. Choices are sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8, sim9."); m->mothurOutEndLine(); abort = true; 
			}
            
            groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "";   }
			else { 
				m->splitAtDash(groups, Groups);	
			}			
			m->setGroups(Groups);
            
            string temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			m->mothurConvert(temp, runs); 

		}

	}
	catch(exception& e) {
		m->errorOut(e, "CooccurrenceCommand", "CooccurrenceCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int CooccurrenceCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		InputData* input = new InputData(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

        ofstream out;
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + "cooccurence.summary";
        m->openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName);  outputTypes["summary"].push_back(outputFileName);
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        out << "metric\tlabel\tScore\tpValue\n";

		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } delete input; out.close(); m->mothurRemove(outputFileName); return 0; }
	
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				getCooccurrence(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
			
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				lookup = input->getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				getCooccurrence(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
			
			if (m->control_pressed) {  outputTypes.clear(); delete input; out.close(); m->mothurRemove(outputFileName); return 0; }

			//get next line to process
			lookup = input->getSharedRAbundVectors();				
		}
		
		if (m->control_pressed) { delete input; out.close(); m->mothurRemove(outputFileName); return 0; }

		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
	
		//run last label if you need to
		if (needToRun == true)  {
			for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }  
			lookup = input->getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			
			getCooccurrence(lookup, out);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
	
        out.close(); 
        
		//reset groups parameter 
		delete input; 

        m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CooccurrenceCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int CooccurrenceCommand::getCooccurrence(vector<SharedRAbundVector*>& thisLookUp, ofstream& out){
	try {
        int numOTUS = thisLookUp[0]->getNumBins();
        vector< vector<int> > initmatrix (thisLookUp.size());
        vector< vector<int> > co_matrix (thisLookUp[0]->getNumBins());
        for (int i = 0; i < thisLookUp[0]->getNumBins(); i++) { co_matrix[i].resize((thisLookUp.size()), 0); }
        for (int i = 0; i < thisLookUp.size(); i++) { initmatrix[i].resize((thisLookUp[i]->getNumBins()), 0); }
        vector<int> columntotal(thisLookUp.size(), 0);
        vector<int> rowtotal(numOTUS, 0);
        
        int rowcount = 0;
        for (int i = 0; i < thisLookUp.size(); i++) {
			for (int j = 0; j < thisLookUp[i]->getNumBins(); j++) {
				if (m->control_pressed) { return 0; }			
				int abund = thisLookUp[i]->getAbundance(j);
				
				if(abund > 0) {
				    initmatrix[i][j] = 1;
                    co_matrix[j][i] = 1;
                    rowcount++;
                    columntotal[j]++;
				}
			}
            rowtotal[i] = rowcount;
            rowcount = 0;
        }
        
        //nrows is ncols of inital matrix. All the functions need this value. They assume the transposition has already taken place and nrows and ncols refer to that matrix.
        //comatrix and initmatrix are still vectors of vectors of ints as in the original script. The abundancevector is only what was read in ie not a co-occurrence matrix!
        int ncols = numOTUS;//rows of inital matrix
        int nrows = thisLookUp.size();//OTUs
        double initscore = 0.0;
        //transpose matrix
        int newmatrows = ncols;
        int newmatcols = nrows;
      
        //swap for transposed matrix
        nrows = newmatrows;//ncols;
        ncols = newmatcols;//nrows;
        
        vector<int> initcolumntotal(ncols, 0);
        vector<int> initrowtotal(nrows, 0);
        vector<double> stats;
               
        TrialSwap2 trial;
        
        initcolumntotal = rowtotal;
        initrowtotal = columntotal;
        trial.update_row_col_totals(co_matrix, rowtotal, columntotal);
        
        if (metric == "cscore")         { initscore = trial.calc_c_score(co_matrix, rowtotal);    }
        else if (metric == "checker")   { initscore = trial.calc_checker(co_matrix, rowtotal);    }
        else if (metric == "vratio")    { initscore = trial.calc_vratio(rowtotal, columntotal);   }
        else if (metric == "combo")     { initscore = trial.calc_combo(co_matrix);                }
        else                            {  m->mothurOut("[ERROR]: No metric selected!\n");  m->control_pressed = true; return 1;            }
        
        m->mothurOut("Initial c score: " + toString(initscore)); m->mothurOutEndLine();
        
        //nullmatrix burn in
        for(int i=0;i<10000;i++) {
            if (m->control_pressed) { return 0; }
            if (matrix == "sim1") {
                trial.sim1(co_matrix);
            }else if (matrix == "sim2") {
                trial.sim2(co_matrix);
            }else if (matrix == "sim3") {
                trial.sim3(initmatrix);
                co_matrix = initmatrix;
            }else if (matrix == "sim4") {
                trial.sim4(columntotal, rowtotal, co_matrix);
            }else if (matrix == "sim5") {
                trial.sim5(initcolumntotal, initrowtotal, initmatrix);
                trial.transpose_matrix(initmatrix,co_matrix);
            }else if (matrix == "sim6") {
                trial.sim6(columntotal, co_matrix);
            }else if (matrix == "sim7") {
                trial.sim7(initcolumntotal, initmatrix);          
                co_matrix = initmatrix;
            }else if (matrix == "sim8") {
                trial.sim8(columntotal, rowtotal, co_matrix);
            }else if (matrix == "sim9") {
                trial.swap_checkerboards (co_matrix);
            }else{
                m->mothurOut("[ERROR]: No model selected! \n");
                m->control_pressed = true;
            }
        }
                
        //run
        for(int i=0;i<runs;i++) {
            if (m->control_pressed) { return 0; }
            //calc metric of nullmatrix
            if (matrix == "sim1") {
                trial.sim1(co_matrix);
            }else if (matrix == "sim2") {
                trial.sim2(co_matrix);
            }else if (matrix == "sim3") {
                trial.sim3(initmatrix);
                co_matrix = initmatrix;
            }else if (matrix == "sim4") {
                trial.sim4(columntotal, rowtotal, co_matrix);
            }else if (matrix == "sim5") {
                trial.sim5(initcolumntotal, initrowtotal, initmatrix);
                trial.transpose_matrix(initmatrix,co_matrix);
            }else if (matrix == "sim6") {
                trial.sim6(columntotal, co_matrix);
            }else if (matrix == "sim7") {
                trial.sim7(initcolumntotal, initmatrix);          
                co_matrix = initmatrix;
            }else if (matrix == "sim8") {
                trial.sim8(columntotal, rowtotal, co_matrix);
            }else if (matrix == "sim9") {
                trial.swap_checkerboards (co_matrix);
            }else{
                 m->mothurOut("[ERROR]: No model selected! \n");
                 m->control_pressed = true;
            }
            //
            //            
            trial.update_row_col_totals(co_matrix, rowtotal, columntotal); 
            
            if (metric == "cscore") { 
                stats.push_back(trial.calc_c_score(co_matrix, rowtotal));
            }else if (metric == "checker") { 
                stats.push_back(trial.calc_checker(co_matrix, rowtotal));
            }else if (metric == "vratio") { 
                stats.push_back(trial.calc_vratio(rowtotal, columntotal));
            }else if (metric == "combo") { 
                stats.push_back(trial.calc_combo(co_matrix));
            }else {
                m->mothurOut("[ERROR]: No metric selected!\n");
                m->control_pressed = true;
                return 1;
            }
            
        }

        double total = 0.0;
        for (int i=0; i<stats.size();i++)   {   total+=stats[i];   }
        
        double nullMean = double (total/(double)stats.size()); 
        
        m->mothurOutEndLine(); m->mothurOut("average metric score: " + toString(nullMean)); m->mothurOutEndLine();
        
        double pvalue = 0.0;
        if (metric == "cscore" || metric == "checker") {    pvalue = trial.calc_pvalue_greaterthan (stats, initscore);   }
        else{   pvalue = trial.calc_pvalue_lessthan (stats, initscore); }
        
        m->mothurOut("pvalue: " + toString(pvalue)); m->mothurOutEndLine();
        out << metric << '\t' << thisLookUp[0]->getLabel() << '\t' << nullMean << '\t' << pvalue << endl;
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CooccurrenceCommand", "Cooccurrence");
		exit(1);
	}
}
//**********************************************************************************************************************


