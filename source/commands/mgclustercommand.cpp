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
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-opti", "opti", "", "", "","",false,false); parameters.push_back(pmethod);
        CommandParameter pinitialize("initialize", "Multiple", "oneotu-singleton", "singleton", "", "", "","",false,false,true); parameters.push_back(pinitialize);
        CommandParameter pmetric("metric", "Multiple", "mcc-sens-spec-tptn-fpfn-tp-tn-fp-fn-f1score-accuracy-ppv-npv-fdr", "mcc", "", "", "","",false,false,true); parameters.push_back(pmetric);
        CommandParameter pmetriccutoff("delta", "Number", "", "0.0001", "", "", "","",false,false,true); parameters.push_back(pmetriccutoff);
        CommandParameter piters("iters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(piters);
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
		helpString += "The mgcluster command parameter options are blast, name, cutoff, precision,   method, metric, initialize, iters, merge, min, length, penalty and adjust. The blast parameter is required.\n";
		helpString += "The mgcluster command reads a blast and name file and clusters the sequences into OPF units similar to the OTUs.\n";
		helpString += "This command outputs a .list, .rabund and .sabund file that can be used with mothur other commands to estimate richness.\n";
		helpString += "The cutoff parameter is used to specify the maximum distance you would like to cluster to. The default is 0.70.\n";
		helpString += "The precision parameter's default value is 100. \n";
		helpString += "The acceptable mgcluster methods are furthest, nearest, average and opti.  If no method is provided then opti is assumed.\n";
		helpString += "The min parameter allows you to specify is you want the minimum or maximum blast score ratio used in calculating the distance. The default is true, meaning you want the minimum.\n";
        helpString += "The iters parameter allow you to set the maxiters for the opticluster method. \n";
        helpString += "The metric parameter allows to select the metric in the opticluster method. Options are Matthews correlation coefficient (mcc), sensitivity (sens), specificity (spec), true positives + true negatives (tptn), false positives + false negatives (fpfn), true positives (tp), true negative (tn), false positive (fp), false negative (fn), f1score (f1score), accuracy (accuracy), positive predictive value (ppv), negative predictive value (npv), false discovery rate (fdr). Default=mcc.\n";
        helpString += "The initialize parameter allows to select the initial randomization for the opticluster method. Options are singleton, meaning each sequence is randomly assigned to its own OTU, or oneotu meaning all sequences are assigned to one otu. Default=singleton.\n";
        helpString += "The delta parameter allows to set the stable value for the metric in the opticluster method (delta=0.0001). \n";
		helpString += "The length parameter is used to specify the minimum overlap required.  The default is 5.\n";
        helpString += "The adjust parameter is used to handle missing distances.  If you set a cutoff, adjust=f by default.  If not, adjust=t by default. Adjust=f, means ignore missing distances and adjust cutoff as needed with the average neighbor method.  Adjust=t, will treat missing distances as 1.0. You can also set the value the missing distances should be set to, adjust=0.5 would give missing distances a value of 0.5.\n";
		helpString += "The penalty parameter is used to adjust the error rate.  The default is 0.10.\n";
		helpString += "The merge parameter allows you to shut off merging based on overlaps and just cluster.  By default merge is true, meaning you want to merge.\n";
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
        else if (type == "steps") {  pattern = "[filename],[clustertag],steps"; }
        else if (type == "sensspec") {  pattern = "[filename],[clustertag],sensspec"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
        outputTypes["steps"] = tempOutNames;
        outputTypes["sensspec"] = tempOutNames;
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
            outputTypes["steps"] = tempOutNames;
            outputTypes["sensspec"] = tempOutNames;
			
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
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "opti"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "opti")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average or opti."); m->mothurOutEndLine(); abort = true; }
            
            metric = validParameter.validFile(parameters, "metric", false);		if (metric == "not found") { metric = "mcc"; }
            
            if ((metric == "mcc") || (metric == "sens") || (metric == "spec") || (metric == "tptn") || (metric == "tp") || (metric == "tn") || (metric == "fp") || (metric == "fn") || (metric == "f1score") || (metric == "accuracy") || (metric == "ppv") || (metric == "npv") || (metric == "fdr") || (metric == "fpfn") ){ }
            else { m->mothurOut("[ERROR]: Not a valid metric.  Valid metrics are mcc, sens, spec, tp, tn, fp, fn, tptn, fpfn, f1score, accuracy, ppv, npv, fdr."); m->mothurOutEndLine(); abort = true; }
            
            initialize = validParameter.validFile(parameters, "initialize", false);		if (initialize == "not found") { initialize = "singleton"; }
            
            if ((initialize == "singleton") || (initialize == "oneotu")){ }
            else { m->mothurOut("[ERROR]: Not a valid initialization.  Valid initializations are singleton and oneotu."); m->mothurOutEndLine(); abort = true; }
            
            temp = validParameter.validFile(parameters, "delta", false);		if (temp == "not found")  { temp = "0.0001"; }
            m->mothurConvert(temp, stableMetric);
            
            temp = validParameter.validFile(parameters, "iters", false);		if (temp == "not found")  { temp = "100"; }
            m->mothurConvert(temp, maxIters);

			temp = validParameter.validFile(parameters, "length", false);			if (temp == "not found") { temp = "5"; }
			m->mothurConvert(temp, length); 
			
			temp = validParameter.validFile(parameters, "penalty", false);			if (temp == "not found") { temp = "0.10"; }
			m->mothurConvert(temp, penalty); 
			
			temp = validParameter.validFile(parameters, "min", false);				if (temp == "not found") { temp = "true"; }
			minWanted = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "merge", false);			if (temp == "not found") { temp = "true"; }
			merge = m->isTrue(temp);
            
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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        fileroot = outputDir + m->getRootName(m->getSimpleName(blastfile));
        tag = "";
        if (method == "furthest")		{ tag = "fn";  }
        else if (method == "nearest")	{ tag = "nn";  }
        else if (method == "average")	{ tag = "an";  }
        else if (method == "opti")      { tag = "opti"; }
        
        if (method == "opti") {  runOptiCluster(); }
        else { runMothurCluster();  }
				
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(listFileName); m->mothurOutEndLine();	outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
		if (countfile == "") {
            m->mothurOut(rabundFileName); m->mothurOutEndLine();	outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            m->mothurOut(sabundFileName); m->mothurOutEndLine();	outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
        }
		m->mothurOutEndLine();
		
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
        }else { mergedList->print(listFile, true); }
        
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
int MGClusterCommand::runOptiCluster(){
    try {
        if (!cutoffSet) {  m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); cutoff = 0.03;  }
        
        string nameOrCount = "";
        string thisNamefile = "";
        map<string, int> counts;
        if (countfile != "") { nameOrCount = "count"; thisNamefile = countfile; CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap(); }
        else if (namefile != "") { nameOrCount = "name"; thisNamefile = namefile; }
        
        string distfile = blastfile;
        
        time_t start = time(NULL);
        
        OptiMatrix matrix;
        matrix.setBlastVariables(length, penalty, minWanted);
        matrix.readFile(distfile, thisNamefile, nameOrCount, "blast", cutoff, false);
        
        ClusterMetric* metricCalc = NULL;
        if (metric == "mcc")             { metricCalc = new MCC();              }
        else if (metric == "sens")       { metricCalc = new Sensitivity();      }
        else if (metric == "spec")       { metricCalc = new Specificity();      }
        else if (metric == "tptn")       { metricCalc = new TPTN();             }
        else if (metric == "tp")         { metricCalc = new TP();               }
        else if (metric == "tn")         { metricCalc = new TN();               }
        else if (metric == "fp")         { metricCalc = new FP();               }
        else if (metric == "fn")         { metricCalc = new FN();               }
        else if (metric == "f1score")    { metricCalc = new F1Score();          }
        else if (metric == "accuracy")   { metricCalc = new Accuracy();         }
        else if (metric == "ppv")        { metricCalc = new PPV();              }
        else if (metric == "npv")        { metricCalc = new NPV();              }
        else if (metric == "fdr")        { metricCalc = new FDR();              }
        else if (metric == "fpfn")       { metricCalc = new FPFN();             }
        
        OptiCluster cluster(&matrix, metricCalc, 0);
        string tag = cluster.getTag();
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        sabundFileName = getOutputFileName("sabund", variables);
        rabundFileName = getOutputFileName("rabund", variables);
        if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        listFileName = getOutputFileName("list", variables);
        string outputName = getOutputFileName("steps", variables);
        outputNames.push_back(outputName); outputTypes["steps"].push_back(outputName);

        m->mothurOutEndLine(); m->mothurOut("Clustering " + distfile); m->mothurOutEndLine();
        
        if (outputDir == "") { outputDir += m->hasPath(distfile); }
        
        ofstream listFile;
        m->openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        
        ofstream outStep;
        m->openOutputFile(outputName, outStep);
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        cluster.initialize(listVectorMetric, true, initialize);
        
        long long numBins = cluster.getNumBins();
        m->mothurOut("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        outStep << "iter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        long long tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        m->mothurOut("0\t0\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(cutoff) + "\t" + toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
        outStep << "0\t0\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t" << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
        for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t"); outStep << results[i] << "\t"; }
        m->mothurOutEndLine();
        outStep << endl;
        
        while ((delta > stableMetric) && (iters < maxIters)) {
            
            int start = time(NULL);
            
            if (m->getControl_pressed()) { break; }
            double oldMetric = listVectorMetric;
            
            cluster.update(listVectorMetric);
            
            delta = abs(oldMetric - listVectorMetric);
            iters++;
            
            results = cluster.getStats(tp, tn, fp, fn);
            numBins = cluster.getNumBins();
            
            m->mothurOut(toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t"+ toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
            outStep << (toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t") << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
            for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t"); outStep << results[i] << "\t"; }
            m->mothurOutEndLine();
            outStep << endl;
        }
        m->mothurOutEndLine(); m->mothurOutEndLine();
        
        list = cluster.getList();
        list->setLabel(toString(cutoff));
        
        if (merge) {
            vector< set<int> > overlap = matrix.getBlastOverlap();
        
            //assign each sequence to bins
            map<string, int> seqToBin;
            for (int i = 0; i < list->getNumBins(); i++) {
                if (m->getControl_pressed()) { break; }
                string bin = list->get(i);
                vector<string> names; m->splitAtComma(bin, names);
                for (int j = 0; j < names.size(); j++) { seqToBin[names[j]] = i; }
            }
            
            //merge overlapping bins
            int mergedBinCount = 0;
            for (int i = 0; i < overlap.size(); i++) {
                set<int> temp = overlap[i]; overlap[i].clear();
                for (set<int>::iterator itOverlap = temp.begin(); itOverlap != temp.end(); itOverlap++) {
                    string firstName = matrix.getOverlapName(i);
                    string secondName = matrix.getOverlapName(*itOverlap);
                    int binKeep = seqToBin[firstName];
                    int binRemove = seqToBin[secondName];
                    
                    if(binKeep != binRemove) {
                        //save names in old bin
                        string bin = list->get(binRemove);
                        
                        //merge bins into name1s bin
                        list->set(binKeep, bin+','+list->get(binKeep));
                        list->set(binRemove, "");
                        
                        vector<string> binNames; m->splitAtComma(bin, binNames);
                        //update binInfo //save name and new bin number
                        for (int k = 0; k < binNames.size(); k++) { seqToBin[binNames[k]] = binKeep; }
                        mergedBinCount++;
                    }
                }
            }
            
            if (mergedBinCount != 0) { m->mothurOut("Merged " + toString(mergedBinCount) + " OTUs based on blast overlap.\n\n"); }
        }
        
        if (!m->getPrintedListHeaders()) { vector<string> temp; m->setListBinLabelsInFile(temp); list->printHeaders(listFile); }
        if(countfile != "") { list->print(listFile, counts); }
        else { list->print(listFile); }
        listFile.close();
        
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	sabundFile);
            m->openOutputFile(rabundFileName,	rabundFile);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
            SAbundVector sabund = list->getSAbundVector();
            sabund.print(sabundFile);
            sabundFile.close();
            
            RAbundVector rabund = list->getRAbundVector();
            rabund.print(rabundFile);
            rabundFile.close();
        }
        delete list;
        
        string sensspecFilename = fileroot+ tag + ".sensspec";
        ofstream sensFile;
        m->openOutputFile(sensspecFilename,	sensFile);
        outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
        
        
        sensFile << "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        
        results = cluster.getStats(tp, tn, fp, fn);
        
        sensFile << cutoff << '\t' << cutoff << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
        for (int i = 0; i < results.size(); i++) {  sensFile << results[i] << '\t'; }
        sensFile << '\n';
        sensFile.close();
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " seconds to cluster."); m->mothurOutEndLine();

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MGClusterCommand", "runOptiCluster");
        exit(1);
    }
}
//**********************************************************************************************************************
int MGClusterCommand::runMothurCluster(){
    try {
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
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        sabundFileName = getOutputFileName("sabund", variables);
        rabundFileName = getOutputFileName("rabund", variables);
        if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        listFileName = getOutputFileName("list", variables);
        
        float previousDist = 0.00000;
        float rndPreviousDist = 0.00000;
        
        time_t start = time(NULL);
        
        //read blastfile - creates sparsematrices for the distances and overlaps as well as a listvector
        //must remember to delete those objects here since readBlast does not
        read = new ReadBlast(blastfile, cutoff, penalty, length, minWanted);
        read->read(nameMap);
        
        list = new ListVector(nameMap->getListVector());
        RAbundVector* rabund = NULL;
        
        if(countfile != "") {
            rabund = new RAbundVector();
            createRabund(ct, list, rabund);
        }else {
            rabund = new RAbundVector(list->getRAbundVector());
        }
        
        if (m->getControl_pressed()) { outputTypes.clear(); delete nameMap; delete read; delete list; delete rabund; return 0; }
        
        
        oldList = *list;
        map<string, int> Seq2Bin;
        map<string, int> oldSeq2Bin;
        
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	sabundFile);
            m->openOutputFile(rabundFileName,	rabundFile);
        }
        m->openOutputFile(listFileName,	listFile);
        list->printHeaders(listFile);
        
        if (m->getControl_pressed()) {
            delete nameMap; delete read; delete list; delete rabund;
            listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();  m->mothurRemove(rabundFileName); m->mothurRemove(sabundFileName); } m->mothurRemove(listFileName);
            outputTypes.clear();
            return 0;
        }
        
        double saveCutoff = cutoff;
        
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
        
        if (m->getControl_pressed()) {
            delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster;
            listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();   m->mothurRemove(rabundFileName); m->mothurRemove(sabundFileName); } m->mothurRemove(listFileName);
            outputTypes.clear();
            return 0;
        }
        
        
        //cluster using cluster classes
        while (distMatrix->getSmallDist() < cutoff && distMatrix->getNNodes() > 0){
            
            if (m->getDebug()) {  cout << "numNodes=" << distMatrix->getNNodes() << " smallDist = " << distMatrix->getSmallDist() << endl; }
            
            cluster->update(cutoff);
            
            if (m->getControl_pressed()) {
                delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster;
                listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();   m->mothurRemove(rabundFileName); m->mothurRemove(sabundFileName); } m->mothurRemove(listFileName);
                outputTypes.clear();
                return 0;
            }
            
            float dist = distMatrix->getSmallDist();
            float rndDist = m->ceilDist(dist, precision);
            
            if(previousDist <= 0.0000 && dist != previousDist){
                oldList.setLabel("unique");
                printData(&oldList, counts);
            }
            else if(rndDist != rndPreviousDist){
                if (merge) {
                    ListVector* temp = mergeOPFs(oldSeq2Bin, rndPreviousDist);
                    
                    if (m->getControl_pressed()) {
                        delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster; delete temp;
                        listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();   m->mothurRemove(rabundFileName); m->mothurRemove(sabundFileName); } m->mothurRemove(listFileName);
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
                
                if (m->getControl_pressed()) {
                    delete nameMap; delete distMatrix; delete list; delete rabund; delete cluster; delete temp;
                    listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();   m->mothurRemove(rabundFileName); m->mothurRemove(sabundFileName); } m->mothurRemove(listFileName);
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
        delete list;
        delete rabund;
        listFile.close();
        
        if (countfile == "") {
            sabundFile.close();
            rabundFile.close();
        }
        if (m->getControl_pressed()) {
            delete nameMap;
            listFile.close(); if (countfile == "") { rabundFile.close(); sabundFile.close();   m->mothurRemove(rabundFileName); m->mothurRemove(sabundFileName); } m->mothurRemove(listFileName);
            outputTypes.clear();
            return 0; 
        }
        
        if (saveCutoff != cutoff) {
            saveCutoff = m->ceilDist(saveCutoff, precision);
            m->mothurOut("changed cutoff to " + toString(cutoff)); m->mothurOutEndLine();
        }
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " seconds to cluster."); m->mothurOutEndLine();
        
        return 0;

    }
    catch(exception& e) {
        m->errorOut(e, "MGClusterCommand", "runMothurCluster");
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
		
		if (overlapMatrix.size() == 0)  {  done = true;  }
		
		while (!done) {
			if (m->getControl_pressed()) {  return newList; }
			
			//get next overlap
			seqDist overlapNode;
			 
            if (count < overlapMatrix.size()) { //do we have another node in the matrix
                overlapNode = overlapMatrix[count];
                count++;
            }else { break; }
			
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
					
                    vector<string> binNames; m->splitAtComma(names, binNames);
					//update binInfo //save name and new bin number
                    for (int i = 0; i < binNames.size(); i++) { binInfo[binNames[i]] = binKeep; }
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


