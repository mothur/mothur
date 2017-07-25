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
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pshared);		
		CommandParameter pmetric("metric", "Multiple", "cscore-checker-combo-vratio", "cscore", "", "", "","",false,false); parameters.push_back(pmetric);
		CommandParameter pmatrix("matrixmodel", "Multiple", "sim1-sim2-sim3-sim4-sim5-sim6-sim7-sim8-sim9", "sim2", "", "", "","",false,false); parameters.push_back(pmatrix);
        CommandParameter pruns("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(pruns);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);

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
		string helpString = "The cooccurrence command calculates four metrics and tests their significance to assess whether presence-absence patterns are different than what one would expect by chance.";
        helpString += "The cooccurrence command parameters are shared, metric, matrixmodel, iters, label and groups.";
        helpString += "The matrixmodel parameter options are sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8 and sim9. Default=sim2";
        helpString += "The metric parameter options are cscore, checker, combo and vratio. Default=cscore";
        helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like analyzed.\n";
        helpString += "The cooccurrence command should be in the following format: \n";
		helpString += "cooccurrence(shared=yourSharedFile) \n";
		helpString += "Example cooccurrence(shared=final.an.shared).\n";
		helpString += "Note: No spaces between parameter labels (i.e. shared), '=' and parameters (i.e.yourShared).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CooccurrenceCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string CooccurrenceCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],cooccurence.summary"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "CooccurrenceCommand", "getOutputPattern");
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
			
			matrix = validParameter.validFile(parameters, "matrixmodel", false);				if (matrix == "not found") { matrix = "sim2"; }
			
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
		
		InputData input(sharedfile, "sharedfile");
		SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

        ofstream out;
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
		string outputFileName = getOutputFileName("summary", variables);
        m->openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName);  outputTypes["summary"].push_back(outputFileName);
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        out << "metric\tlabel\tScore\tzScore\tstandardDeviation\tnp_Pvalue\n";

		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
            if (m->control_pressed) { delete lookup; out.close(); m->mothurRemove(outputFileName); return 0; }
	
			if(allLines == 1 || labels.count(lookup->getLabel()) == 1){

				m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
				
				getCooccurrence(lookup, out);
				
				processedLabels.insert(lookup->getLabel());
				userLabels.erase(lookup->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup->getLabel();
			
				delete lookup;
				lookup = input.getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
				getCooccurrence(lookup, out);
				
				processedLabels.insert(lookup->getLabel());
				userLabels.erase(lookup->getLabel());
				
				//restore real lastlabel to save below
				lookup->setLabels(saveLabel);
			}
			
			lastLabel = lookup->getLabel();
			delete lookup;
			
			if (m->control_pressed) {  outputTypes.clear(); out.close(); m->mothurRemove(outputFileName); return 0; }

			//get next line to process
			lookup = input.getSharedRAbundVectors();
		}
		
		if (m->control_pressed) { out.close(); m->mothurRemove(outputFileName); return 0; }

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
			delete lookup;
			lookup = input.getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
			
			getCooccurrence(lookup, out);
			
			delete lookup;
		}
	
        out.close(); 
        
		//reset groups parameter
        m->clearGroups(); 

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

int CooccurrenceCommand::getCooccurrence(SharedRAbundVectors* thisLookUp, ofstream& out){
    try {
        int numOTUS = thisLookUp->getNumBins();
        
        if(numOTUS < 2) {
            m->mothurOut("Not enough OTUs for co-occurrence analysis, skipping"); m->mothurOutEndLine();
            return 0;
        }
        
        vector< vector<int> > co_matrix; co_matrix.resize(thisLookUp->getNumBins());
        for (int i = 0; i < thisLookUp->getNumBins(); i++) { co_matrix[i].resize((thisLookUp->size()), 0); }
        vector<int> columntotal; columntotal.resize(thisLookUp->size(), 0);
        vector<int> rowtotal; rowtotal.resize(numOTUS, 0);
        
        
        vector<RAbundVector*> temp;
        for (int i = 0; i < thisLookUp->size(); i++) { //nrows in the shared file
            for (int j = 0; j < temp[i]->getNumBins(); j++) { //cols of original shared file
                if (m->control_pressed) { return 0; }
                int abund = temp[i]->get(j);
                
                if(abund > 0) {
                    co_matrix[j][i] = 1;
                    rowtotal[j]++;
                    columntotal[i]++;
                }
            }
        }
        for (int i = 0; i < temp.size(); i++) { delete temp[i]; }
        
        //nrows is ncols of inital matrix. All the functions need this value. They assume the transposition has already taken place and nrows and ncols refer to that matrix.
        //comatrix and initmatrix are still vectors of vectors of ints as in the original script. The abundancevector is only what was read in ie not a co-occurrence matrix!
        int nrows = numOTUS;//rows of inital matrix
        int ncols = thisLookUp->size();//groups
        double initscore = 0.0;
       
        vector<double> stats;
        vector<double> probabilityMatrix; probabilityMatrix.resize(ncols * nrows, 0);
        vector<vector<int> > nullmatrix(nrows, vector<int>(ncols, 0));
       
        TrialSwap2 trial;
        
        int n = accumulate( columntotal.begin(), columntotal.end(), 0 );
        
        //============================================================
        
        //generate a probability matrix. Only do this once.
        float start = 0.0;
        
        if (matrix == "sim1") {
            for(int i=0;i<nrows;i++) {
                for(int j=0;j<ncols;j++) {
                    probabilityMatrix[ncols * i + j] = start + 1/double(nrows*ncols);
                    start = start + 1/double(nrows*ncols);
                }
            }
        }
        //don't need a prob matrix because we just shuffle the rows, may use this in the future
        else if (matrix == "sim2") { }
//            for(int i=0;i<nrows;i++) {
//                start = 0.0;
//                for(int j=0;j<ncols;j++) {
//                    probabilityMatrix[ncols * i + j] = start + 1/double(ncols);
//                    start = start + 1/double(ncols);
//                }
//            }
//        }
        
        else if (matrix == "sim3") {
            for(int j=0;j<ncols;j++) {
                start = 0.0;
                for(int i=0;i<nrows;i++) {
                    probabilityMatrix[ncols * i + j] = start + 1/double(nrows);
                    start = start + 1/double(nrows);
                }
            }
        }
        
        else if (matrix == "sim4") {
            for(int i=0;i<nrows;i++) {
                start = 0.0;
                for(int j=0;j<ncols;j++) {
                    probabilityMatrix[ncols * i + j] = start + columntotal[j]/double(n);
                    start = start + columntotal[j]/double(n);
                }
            }
        }
        
        else if (matrix == "sim5") {
            for(int j=0;j<ncols;j++) {
                start = 0.0;
                for(int i=0;i<nrows;i++) {
                    probabilityMatrix[ncols * i + j] = start + rowtotal[i]/double(n);
                    start = start + rowtotal[i]/double(n);
                }
            }
        }
        
        else if (matrix == "sim6") {
            for(int i=0;i<nrows;i++) {
                for(int j=0;j<ncols;j++) {
                    probabilityMatrix[ncols * i + j] = start + columntotal[j]/double(n*nrows);
                    start = start + columntotal[j]/double(n*nrows);
                }
            }
        }
        
        
        else if (matrix == "sim7") {
            for(int i=0;i<nrows;i++) {
                for(int j=0;j<ncols;j++) {
                    probabilityMatrix[ncols * i + j] = start + rowtotal[i]/double(n*ncols);
                    start = start + rowtotal[i]/double(n*ncols);
                }
            }
        }
        
        else if (matrix == "sim8") {
            for(int i=0;i<nrows;i++) {
                for(int j=0;j<ncols;j++) {
                    probabilityMatrix[ncols * i + j] = start + (rowtotal[i]*columntotal[j])/double(n*n);
                    start = start + (rowtotal[i]*columntotal[j])/double(n*n);
                }
            }
        }
        else if (matrix == "sim9" || matrix == "sim2") { }
        else {
            m->mothurOut("[ERROR]: No model selected! \n");
            m->control_pressed = true;
        }
        
        
        //co_matrix is the transposed shared file, initmatrix is the original shared file
        if (metric == "cscore") { initscore = trial.calc_c_score(co_matrix, rowtotal, ncols, nrows); }
        else if (metric == "checker") { initscore = trial.calc_checker(co_matrix, rowtotal, ncols, nrows); }
        else if (metric == "vratio") { initscore = trial.calc_vratio(nrows, ncols, rowtotal, columntotal); }
        else if (metric == "combo") { initscore = trial.calc_combo(nrows, ncols, co_matrix); }
        else { m->mothurOut("[ERROR]: No metric selected!\n"); m->control_pressed = true; return 1; }
        
        m->mothurOut("Initial c score: " + toString(initscore)); m->mothurOutEndLine();
        
        double previous;
        double current;
        double randnum;
        int count;

        //burn-in for sim9    
        if(matrix == "sim9") {
            for(int i=0;i<10000;i++) trial.swap_checkerboards (co_matrix, ncols, nrows);
        }

        //populate null matrix from probability matrix, do this a lot.
        for(int k=0;k<runs;k++){
            nullmatrix.clear();
            //zero-fill the null matrix
            nullmatrix.assign(nrows, vector<int>(ncols, 0));
            
            if(matrix == "sim1" || matrix == "sim6" || matrix == "sim8" || matrix == "sim7") {
                count = 0;
                while(count < n) {
                    if (m->control_pressed) { return 0; }
                nextnum2:
                    previous = 0.0;
                    randnum = m->getRandomDouble0to1();
                    for(int i=0;i<nrows;i++) {
                        for(int j=0;j<ncols;j++) {
                            current = probabilityMatrix[ncols * i + j];
                            if(randnum <= current && randnum > previous) {
                                nullmatrix[i][j] = 1;
                                count++;
                                if (count > n) break;
                                else
                                    goto nextnum2;
                            }
                            previous = current;
                        }
                    }
                }
            }
            
            else if (matrix == "sim2") {
                for(int i=0;i<nrows;i++) {
                    m->mothurRandomShuffle(co_matrix[i]);
                }
                //do this for the scoring since those all have nullmatrix as a parameter
                //nullmatrix gets cleared at the begining of each run
                nullmatrix = co_matrix;
            }
            
            else if(matrix == "sim4") {
                for(int i=0;i<nrows;i++) {
                    count = 0;
                    while(count < rowtotal[i]) {
                        previous = 0.0;
                        if (m->control_pressed) { return 0; }
                        randnum = m->getRandomDouble0to1();
                        for(int j=0;j<ncols;j++) {
                            current = probabilityMatrix[ncols * i + j];
                            if(randnum <= current && randnum > previous && nullmatrix[i][j] != 1) {
                                nullmatrix[i][j] = 1;
                                count++;
                                previous = 0.0;
                                break;
                            }
                            previous = current;
                        }
                    }
                }
            }
            
            else if(matrix == "sim3" || matrix == "sim5") {
                //columns
                for(int j=0;j<ncols;j++) {
                    count = 0;
                    while(count < columntotal[j]) {
                        if (m->control_pressed) { return 0; }
                        randnum = m->getRandomDouble0to1();
                        for(int i=0;i<nrows;i++) {
                            current = probabilityMatrix[ncols * i + j];
                            if(randnum <= current && randnum > previous && nullmatrix[i][j] != 1) {
                                nullmatrix[i][j] = 1;
                                count++;
                                previous = 0.0;
                                break;
                            }
                            previous = current;
                        }
                    }
                }
            }
            
            //swap_checkerboards takes the original matrix and swaps checkerboards
            else if(matrix == "sim9") {
                trial.swap_checkerboards (co_matrix, ncols, nrows);
                nullmatrix = co_matrix;
            }
            else {
                m->mothurOut("[ERROR]: No null model selected!\n\n"); m->control_pressed = true;
                return 1;
            }
            
            //run metric on null matrix and add score to the stats vector
            if (metric == "cscore"){
                stats.push_back(trial.calc_c_score(nullmatrix, rowtotal, ncols, nrows));
            }
            else if (metric == "checker") {
                stats.push_back(trial.calc_checker(nullmatrix, rowtotal, ncols, nrows));
            }
            else if (metric == "vratio") {
                stats.push_back(trial.calc_vratio(nrows, ncols, rowtotal, columntotal));
            }
            else if (metric == "combo") {
                stats.push_back(trial.calc_combo(nrows, ncols, nullmatrix));
            }
            else {
                m->mothurOut("[ERROR]: No metric selected!\n\n"); m->control_pressed = true;
                return 1;
            }
            
        }
        
        
        
        double total = 0.0;
        for (int i=0; i<stats.size();i++) { total+=stats[i]; }
        
        double nullMean = double (total/(double)stats.size());
        
        m->mothurOutEndLine(); m->mothurOut("average metric score: " + toString(nullMean)); m->mothurOutEndLine();
        
        //calc_p_value is not a statistical p-value, it's just the average that are either > or < the initscore.
        //All it does is show what is expected in a competitively structured community
        //zscore is output so p-value can be looked up in a ztable
        double pvalue = 0.0;
        if (metric == "cscore" || metric == "checker") { pvalue = trial.calc_pvalue_greaterthan (stats, initscore); }
        else{ pvalue = trial.calc_pvalue_lessthan (stats, initscore); }

        double sd = trial.getSD(runs, stats, nullMean);

        double zscore = trial.get_zscore(sd, nullMean, initscore);
        
        m->mothurOut("zscore: " + toString(zscore)); m->mothurOutEndLine();
        m->mothurOut("standard deviation: " + toString(sd)); m->mothurOutEndLine();
        m->mothurOut("non-parametric p-value: " + toString(pvalue)); m->mothurOutEndLine();
        out << metric << '\t' << thisLookUp->getLabel() << '\t' << nullMean << '\t' << zscore << '\t' << sd << '\t' << pvalue << endl;
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CooccurrenceCommand", "Cooccurrence");
        exit(1);
    }
}
//**********************************************************************************************************************


