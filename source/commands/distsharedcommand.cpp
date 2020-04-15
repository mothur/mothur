/*
 *  distsharedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "distsharedcommand.h"
#include "subsample.h"

//**********************************************************************************************************************
vector<string> DistSharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","phylip",false,true,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter psubsample("subsample", "String", "", "", "", "", "","",false,false); parameters.push_back(psubsample);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pcalc("calc", "Multiple", "sharedsobs-sharedchao-sharedace-jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-kstest-sharednseqs-ochiai-anderberg-kulczynski-kulczynskicody-lennon-morisitahorn-braycurtis-whittaker-odum-canberra-structeuclidean-structchord-hellinger-manhattan-structpearson-soergel-spearman-structkulczynski-speciesprofile-hamming-structchi2-gower-memchi2-memchord-memeuclidean-mempearson-jsd-rjsd", "jclass-thetayc", "", "", "","",true,false,true); parameters.push_back(pcalc);
		CommandParameter poutput("output", "Multiple", "lt-square-column", "lt", "", "", "","",false,false); parameters.push_back(poutput);
        CommandParameter pmode("mode", "Multiple", "average-median", "average", "", "", "","",false,false); parameters.push_back(pmode);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
        CommandParameter pwithreplacement("withreplacement", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pwithreplacement);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> tempOutNames;
        outputTypes["phylip"] = tempOutNames;
        
        abort = false; calledHelp = false;   allLines = true;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DistSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DistSharedCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The dist.shared command parameters are shared, groups, calc, output, processors, subsample, iters, mode, and label.  shared is a required, unless you have a valid current file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included used.\n";
		helpString += "The group names are separated by dashes. The label parameter allows you to select what distance levels you would like distance matrices created for, and is also separated by dashes.\n";
        helpString += "The iters parameter allows you to choose the number of times you would like to run the subsample.\n";
        helpString += "The subsample parameter allows you to enter the size pergroup of the sample or you can set subsample=T and mothur will use the size of your smallest group.\n";
        helpString += "The withreplacement parameter allows you to indicate you want to subsample your data allowing for the same read to be included multiple times. Default=f. \n";
		helpString += "The dist.shared command should be in the following format: dist.shared(groups=yourGroups, calc=yourCalcs, label=yourLabels).\n";
		helpString += "The output parameter allows you to specify format of your distance matrix. Options are lt, column and square. The default is lt.\n";
        helpString += "The mode parameter allows you to specify if you want the average or the median values reported when subsampling. Options are average, and median. The default is average.\n";
		helpString += "Example dist.shared(groups=A-B-C, calc=jabund-sorabund).\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The default value for calc is jclass and thetayc.\n";
		helpString += validCalculator.printCalc("matrix");
		helpString += "The dist.shared command outputs a .dist file for each calculator you specify at each distance you choose.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DistSharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string DistSharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "phylip") {  pattern = "[filename],[calc],[distance],[outputtag],dist-[filename],[calc],[distance],[outputtag],[tag2],dist"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "DistSharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
DistSharedCommand::DistSharedCommand(string option)  {
	try {

		if(option == "help") {  help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters  = parser.getParameters();
			
			ValidParameters validParameter;
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not found") { 			
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required.\n");  abort = true; }
			}else if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else { current->setSharedFile(sharedfile); }
			
			 
            if (outputdir == ""){ outputdir += util.hasPath(sharedfile); }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			output = validParameter.valid(parameters, "output");		if(output == "not found"){	output = "lt"; }
			if ((output != "lt") && (output != "square") && (output != "column")) { m->mothurOut(output + " is not a valid output form. Options are lt, column and square. I will use lt.\n"); output = "lt"; }
            
            mode = validParameter.valid(parameters, "mode");		if(mode == "not found"){	mode = "average"; }
			if ((mode != "average") && (mode != "median")) { m->mothurOut(mode + " is not a valid mode. Options are average and medina. I will use average.\n");  output = "average"; }
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
            }
			
			string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
				
			calc = validParameter.valid(parameters, "calc");
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			util.splitAtDash(calc, Estimators);
			if (util.inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}
            
            temp = validParameter.valid(parameters, "iters");			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
            
            temp = validParameter.valid(parameters, "subsample");		if (temp == "not found") { temp = "F"; }
			if (util.isNumeric1(temp)) { util.mothurConvert(temp, subsampleSize); subsample = true; }
            else {  
                if (util.isTrue(temp)) { subsample = true; subsampleSize = -1; }  //we will set it to smallest group later 
                else { subsample = false; }
            }
            
            if (subsample == false) { iters = 1; }
            
            temp = validParameter.valid(parameters, "withreplacement");		if (temp == "not found"){	temp = "f";		}
            withReplacement = util.isTrue(temp);
            
        }
	}
	catch(exception& e) {
		m->errorOut(e, "DistSharedCommand", "DistSharedCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

DistSharedCommand::~DistSharedCommand(){}

//**********************************************************************************************************************

int DistSharedCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
	
        time_t start = time(NULL);
        
		InputData input(sharedfile, "sharedfile", Groups);
		set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
					
        if (lookup->size() < 2) { m->mothurOut("[ERROR]: You have not provided enough valid groups.  I cannot run the command.\n");  delete lookup; return 0;}
        
        if (subsample) { 
            if (subsampleSize == -1) { //user has not set size, set size = smallest samples size
                subsampleSize = lookup->getNumSeqsSmallestGroup();
                m->mothurOut("\nSetting sample size to " + toString(subsampleSize) + ".\n\n");
            }else {
                lookup->removeGroups(subsampleSize);
                Groups = lookup->getNamesGroups();
            }
            
            if (lookup->size() < 2) { m->mothurOut("[ERROR]: You have not provided enough valid groups.  I cannot run the command.\n"); m->setControl_pressed(true);  return 0; }
        }
		numGroups = lookup->size();
        
        if (m->getControl_pressed()) { delete lookup;  return 0;  }
        
        while (lookup != NULL) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            createProcesses(lookup); delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
				
		if (m->getControl_pressed()) { outputTypes.clear();  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  return 0;  }
		
		//set phylip file as new current phylipfile
		string currentName = "";
		itTypes = outputTypes.find("phylip");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setPhylipFile(currentName);  }
		}
		
        m->mothurOut("\nIt took " + toString(time(NULL) - start) + " seconds to run dist.shared.\n");
        
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "DistSharedCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
void DistSharedCommand::printDists(ostream& out, vector< vector<double> >& simMatrix, vector<string> groupNames) {
    try {
        
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        if (output == "lt") {
            out << simMatrix.size() << endl;
            for (int b = 0; b < simMatrix.size(); b++)	{
                out << groupNames[b];
                for (int n = 0; n < b; n++)	{
                    out  << '\t' << simMatrix[b][n];
                }
                out << endl;
            }
        }else if (output == "column") {
            for (int b = 0; b < simMatrix.size(); b++)	{
                for (int n = 0; n < b; n++)	{
                    out << groupNames[b] << '\t' << groupNames[n] << '\t' << simMatrix[b][n] << endl;
                }
            }
        }else{
            out << simMatrix.size() << endl;
            for (int b = 0; b < simMatrix.size(); b++)	{
                out << groupNames[b];
                for (int n = 0; n < simMatrix[b].size(); n++)	{
                    out << '\t' << simMatrix[b][n];
                }
                out << endl;
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "DistSharedCommand", "printSims");
        exit(1);
    }
}
/**************************************************************************************************/
int driver(vector<SharedRAbundVector*>& thisLookup, vector< vector<seqDist> >& calcDists, vector<Calculator*> matrixCalculators, MothurOut* m) {
    try {
        vector<SharedRAbundVector*> subset;
        
        for (int k = 0; k < thisLookup.size(); k++) { // pass cdd each set of groups to compare
            
            for (int l = 0; l < k; l++) {
                
                if (k != l) { //we dont need to similarity of a groups to itself
                    subset.clear(); //clear out old pair of sharedrabunds
                    //add new pair of sharedrabunds
                    subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]);
                    
                    for(int i=0;i<matrixCalculators.size();i++) {
                        
                        //if this calc needs all groups to calculate the pair load all groups
                        if (matrixCalculators[i]->getNeedsAll()) {
                            //load subset with rest of lookup for those calcs that need everyone to calc for a pair
                            for (int w = 0; w < thisLookup.size(); w++) {
                                if ((w != k) && (w != l)) { subset.push_back(thisLookup[w]); }
                            }
                        }
                        
                        vector<double> tempdata = matrixCalculators[i]->getValues(subset); //saves the calculator outputs
                        
                        if (m->getControl_pressed()) { return 1; }
                        
                        seqDist temp(l, k, tempdata[0]);
                        calcDists[i].push_back(temp);
                    }
                }
            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "DistSharedCommand", "driver");
        exit(1);
    }
}

/***********************************************************/
int process(distSharedData* params){
	try {
        vector<Calculator*> matrixCalculators;
        ValidCalculators validCalculator;
        for (int i=0; i<params->Estimators.size(); i++) {
            if (validCalculator.isValidCalculator("matrix", params->Estimators[i]) ) {
                if (params->Estimators[i] == "sharedsobs") {
                    matrixCalculators.push_back(new SharedSobsCS());
                }else if (params->Estimators[i] == "sharedchao") {
                    matrixCalculators.push_back(new SharedChao1());
                }else if (params->Estimators[i] == "sharedace") {
                    matrixCalculators.push_back(new SharedAce());
                }else if (params->Estimators[i] == "jabund") {
                    matrixCalculators.push_back(new JAbund());
                }else if (params->Estimators[i] == "sorabund") {
                    matrixCalculators.push_back(new SorAbund());
                }else if (params->Estimators[i] == "jclass") {
                    matrixCalculators.push_back(new Jclass());
                }else if (params->Estimators[i] == "sorclass") {
                    matrixCalculators.push_back(new SorClass());
                }else if (params->Estimators[i] == "jest") {
                    matrixCalculators.push_back(new Jest());
                }else if (params->Estimators[i] == "sorest") {
                    matrixCalculators.push_back(new SorEst());
                }else if (params->Estimators[i] == "thetayc") {
                    matrixCalculators.push_back(new ThetaYC());
                }else if (params->Estimators[i] == "thetan") {
                    matrixCalculators.push_back(new ThetaN());
                }else if (params->Estimators[i] == "kstest") {
                    matrixCalculators.push_back(new KSTest());
                }else if (params->Estimators[i] == "sharednseqs") {
                    matrixCalculators.push_back(new SharedNSeqs());
                }else if (params->Estimators[i] == "ochiai") {
                    matrixCalculators.push_back(new Ochiai());
                }else if (params->Estimators[i] == "anderberg") {
                    matrixCalculators.push_back(new Anderberg());
                }else if (params->Estimators[i] == "kulczynski") {
                    matrixCalculators.push_back(new Kulczynski());
                }else if (params->Estimators[i] == "kulczynskicody") {
                    matrixCalculators.push_back(new KulczynskiCody());
                }else if (params->Estimators[i] == "lennon") {
                    matrixCalculators.push_back(new Lennon());
                }else if (params->Estimators[i] == "morisitahorn") {
                    matrixCalculators.push_back(new MorHorn());
                }else if (params->Estimators[i] == "braycurtis") {
                    matrixCalculators.push_back(new BrayCurtis());
                }else if (params->Estimators[i] == "whittaker") {
                    matrixCalculators.push_back(new Whittaker());
                }else if (params->Estimators[i] == "odum") {
                    matrixCalculators.push_back(new Odum());
                }else if (params->Estimators[i] == "canberra") {
                    matrixCalculators.push_back(new Canberra());
                }else if (params->Estimators[i] == "structeuclidean") {
                    matrixCalculators.push_back(new StructEuclidean());
                }else if (params->Estimators[i] == "structchord") {
                    matrixCalculators.push_back(new StructChord());
                }else if (params->Estimators[i] == "hellinger") {
                    matrixCalculators.push_back(new Hellinger());
                }else if (params->Estimators[i] == "manhattan") {
                    matrixCalculators.push_back(new Manhattan());
                }else if (params->Estimators[i] == "structpearson") {
                    matrixCalculators.push_back(new StructPearson());
                }else if (params->Estimators[i] == "soergel") {
                    matrixCalculators.push_back(new Soergel());
                }else if (params->Estimators[i] == "spearman") {
                    matrixCalculators.push_back(new Spearman());
                }else if (params->Estimators[i] == "structkulczynski") {
                    matrixCalculators.push_back(new StructKulczynski());
                }else if (params->Estimators[i] == "speciesprofile") {
                    matrixCalculators.push_back(new SpeciesProfile());
                }else if (params->Estimators[i] == "hamming") {
                    matrixCalculators.push_back(new Hamming());
                }else if (params->Estimators[i] == "structchi2") {
                    matrixCalculators.push_back(new StructChi2());
                }else if (params->Estimators[i] == "gower") {
                    matrixCalculators.push_back(new Gower());
                }else if (params->Estimators[i] == "memchi2") {
                    matrixCalculators.push_back(new MemChi2());
                }else if (params->Estimators[i] == "memchord") {
                    matrixCalculators.push_back(new MemChord());
                }else if (params->Estimators[i] == "memeuclidean") {
                    matrixCalculators.push_back(new MemEuclidean());
                }else if (params->Estimators[i] == "mempearson") {
                    matrixCalculators.push_back(new MemPearson());
                }else if (params->Estimators[i] == "jsd") {
                    matrixCalculators.push_back(new JSD());
                }else if (params->Estimators[i] == "rjsd") {
                    matrixCalculators.push_back(new RJSD());
                }
            }
        }
        
        //if the users entered no valid calculators don't execute command
        if (matrixCalculators.size() == 0) { params->m->mothurOut("No valid calculators.\n");  return 0; }
        params->Estimators.clear();
        for (int i=0; i<matrixCalculators.size(); i++) { params->Estimators.push_back(matrixCalculators[i]->getName()); }
        
        vector< vector<seqDist>  > calcDists; calcDists.resize(matrixCalculators.size()); 		
        SubSample sample;
        for (int thisIter = 0; thisIter < params->numIters; thisIter++) {
            SharedRAbundVectors* thisItersLookup = new SharedRAbundVectors(*params->thisLookup);
            vector<string> namesOfGroups = thisItersLookup->getNamesGroups();
            
            time_t start = time(NULL);
            
            if (params->subsample) {
                if (params->withReplacement)    {  sample.getSampleWithReplacement(thisItersLookup, params->subsampleSize);     }
                else                            {  sample.getSample(thisItersLookup, params->subsampleSize);                    }
            }
            if (params->m->getDebug()) { params->m->mothurOut("\nIt took " + toString(time(NULL) - start) + " seconds to subsample the shared file.\n");  }
            
            //params->m->mothurOut(toString(thisIter) + " It took " + toString(time(NULL) - start) + " seconds to subsample the shared file.\n");
            
            vector<SharedRAbundVector*> thisItersRabunds = thisItersLookup->getSharedRAbundVectors();
            vector<string> thisItersGroupNames = params->thisLookup->getNamesGroups();
            
            start = time(NULL);
            driver(thisItersRabunds, calcDists, matrixCalculators, params->m);
            if (params->m->getDebug()) { params->m->mothurOut("\nIt took " + toString(time(NULL) - start) + " seconds to calc dist for shared file.\n");  }
            
            //params->m->mothurOut(toString(thisIter) + " It took " + toString(time(NULL) - start) + " seconds to calc dist for shared file.\n");
            
            for (int i = 0; i < thisItersRabunds.size(); i++) { delete thisItersRabunds[i]; }
            if (params->subsample){
                if((thisIter+1) % 100 == 0){	params->m->mothurOutJustToScreen(toString(thisIter+1)+"\n"); 		}
                params->calcDistsTotals.push_back(calcDists);
                for (int i = 0; i < calcDists.size(); i++) {
                    for (int j = 0; j < calcDists[i].size(); j++) {
                        if (params->m->getDebug()) {  params->m->mothurOut("[DEBUG]: Results: iter = " + toString(thisIter) + ", " + thisItersGroupNames[calcDists[i][j].seq1] + " - " + thisItersGroupNames[calcDists[i][j].seq2] + " distance = " + toString(calcDists[i][j].dist) + ".\n");  }
                    } 
                }
            }else { //print results for whole dataset
                for (int i = 0; i < calcDists.size(); i++) {
                    if (params->m->getControl_pressed()) { break; }
                    
                    //initialize matrix
                    vector< vector<double> > matrix; //square matrix to represent the distance
                    matrix.resize(thisItersLookup->size());
                    for (int k = 0; k < thisItersLookup->size(); k++) {  matrix[k].resize(thisItersLookup->size(), 0.0); }
                    
                    for (int j = 0; j < calcDists[i].size(); j++) {
                        int row = calcDists[i][j].seq1;
                        int column = calcDists[i][j].seq2;
                        double dist = calcDists[i][j].dist;
                        
                        matrix[row][column] = dist;
                        matrix[column][row] = dist;
                    }
                    params->matrices.push_back(matrix);
                }
            }
            for (int i = 0; i < calcDists.size(); i++) {  calcDists[i].clear(); }
            delete thisItersLookup;
		}
		if((params->numIters) % 100 != 0){	params->m->mothurOutJustToScreen(toString(params->numIters)+"\n"); 		}
        
		return 0;
	}
	catch(exception& e) {
		params->m->errorOut(e, "DistSharedCommand", "process");
		exit(1);
	}
}
/***********************************************************/
int DistSharedCommand::createProcesses(SharedRAbundVectors*& thisLookup){
    try {
        
        vector<string> groupNames = thisLookup->getNamesGroups();
        
        vector<int> lines;
        if (processors > (iters)) { processors = iters; }
        
        //figure out how many sequences you have to process
        int numItersPerProcessor = (iters) / processors;
        for (int i = 0; i < processors; i++) {
            if(i == (processors - 1)){	numItersPerProcessor = (iters) - i * numItersPerProcessor; 	}
            lines.push_back(numItersPerProcessor);
        }
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<distSharedData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            
            //make copy of lookup so we don't get access violations
            SharedRAbundVectors* newLookup = new SharedRAbundVectors(*thisLookup);
            distSharedData* dataBundle = new distSharedData(lines[i+1], false, subsample, subsampleSize, withReplacement, Estimators, newLookup);
            
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(process, dataBundle));
        }
        
        //make copy of lookup so we don't get access violations
        SharedRAbundVectors* newLookup = new SharedRAbundVectors(*thisLookup);
        distSharedData* dataBundle = new distSharedData(lines[0], true, subsample, subsampleSize, withReplacement, Estimators, newLookup);
        process(dataBundle);
        delete newLookup;
        
        Estimators.clear(); Estimators = dataBundle->Estimators;
        
        if (!subsample) {
            map<string, string> variables;
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
            variables["[distance]"] = thisLookup->getLabel();
            variables["[tag2]"] = "";
            variables["[outputtag]"] = output;
            
            /// fix to print out matrices for each calc - only main does this
            for (int i = 0; i < Estimators.size(); i++) {
                variables["[calc]"] = Estimators[i];
                string distFileName = getOutputFileName("phylip",variables);
                outputNames.push_back(distFileName); outputTypes["phylip"].push_back(distFileName);
                
                ofstream outDist; util.openOutputFile(distFileName, outDist);
                outDist.setf(ios::fixed, ios::floatfield); outDist.setf(ios::showpoint);
                
                printDists(outDist, dataBundle->matrices[i], groupNames); outDist.close();
            }
        }
        vector< vector< vector<seqDist> > > calcDistsTotals = dataBundle->calcDistsTotals;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            //get calcDistsTotal info - one entry per iter
            for (int j = 0; j < data[i]->calcDistsTotals.size(); j++) { calcDistsTotals.push_back(data[i]->calcDistsTotals[j]); }
            
            delete data[i]->thisLookup;
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
    
        if (subsample) {
            //we need to find the average distance and standard deviation for each groups distance
            vector< vector<seqDist>  > calcAverages = util.getAverages(calcDistsTotals, mode);
            
            //find standard deviation
            vector< vector<seqDist>  > stdDev = util.getStandardDeviation(calcDistsTotals, calcAverages);
            
            //print results
            for (int i = 0; i < Estimators.size(); i++) {
                vector< vector<double> > matrix; //square matrix to represent the distance
                matrix.resize(thisLookup->size());
                for (int k = 0; k < thisLookup->size(); k++) {  matrix[k].resize(thisLookup->size(), 0.0); }
                
                vector< vector<double> > stdmatrix; //square matrix to represent the stdDev
                stdmatrix.resize(thisLookup->size());
                for (int k = 0; k < thisLookup->size(); k++) {  stdmatrix[k].resize(thisLookup->size(), 0.0); }
                
                
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    int row = calcAverages[i][j].seq1;
                    int column = calcAverages[i][j].seq2;
                    float dist = calcAverages[i][j].dist;
                    float stdDist = stdDev[i][j].dist;
                    
                    matrix[row][column] = dist;
                    matrix[column][row] = dist;
                    stdmatrix[row][column] = stdDist;
                    stdmatrix[column][row] = stdDist;
                }
                
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
                variables["[distance]"] = thisLookup->getLabel();
                variables["[outputtag]"] = output;
                variables["[tag2]"] = "ave";
                variables["[calc]"] = Estimators[i];
                string distFileName = getOutputFileName("phylip",variables);
                outputNames.push_back(distFileName); outputTypes["phylip"].push_back(distFileName);
                //set current phylip file to average distance matrix
                current->setPhylipFile(distFileName);
                ofstream outAve;
                util.openOutputFile(distFileName, outAve);
                outAve.setf(ios::fixed, ios::floatfield); outAve.setf(ios::showpoint);
                
                printDists(outAve, matrix, groupNames);
                
                outAve.close();
                
                variables["[tag2]"] = "std";
                distFileName = getOutputFileName("phylip",variables);
                outputNames.push_back(distFileName); outputTypes["phylip"].push_back(distFileName);
                ofstream outSTD;
                util.openOutputFile(distFileName, outSTD);
                outSTD.setf(ios::fixed, ios::floatfield); outSTD.setf(ios::showpoint);
                
                printDists(outSTD, stdmatrix, thisLookup->getNamesGroups());
                
                outSTD.close();
            }
        }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "DistSharedCommand", "createProcesses");
        exit(1);
    }
}
/***********************************************************/



