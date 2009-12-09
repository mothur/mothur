/*
 *  classifyseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/2/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "classifyseqscommand.h"
#include "sequence.hpp"
#include "bayesian.h"
#include "phylotree.h"
#include "knn.h"

//**********************************************************************************************************************

ClassifySeqsCommand::ClassifySeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"template","fasta","search","ksize","method","processors","taxonomy","match","mismatch","gapopen","gapextend","numwanted","cutoff","probs"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			templateFileName = validParameter.validFile(parameters, "template", true);
			if (templateFileName == "not found") { 
				mothurOut("template is a required parameter for the classify.seqs command."); 
				mothurOutEndLine();
				abort = true; 
			}
			else if (templateFileName == "not open") { abort = true; }	
			
			fastaFileName = validParameter.validFile(parameters, "fasta", true);
			if (fastaFileName == "not found") { 
				mothurOut("fasta is a required parameter for the classify.seqs command."); 
				mothurOutEndLine();
				abort = true; 
			}
			else if (fastaFileName == "not open") { abort = true; }	
			
			taxonomyFileName = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyFileName == "not found") { 
				mothurOut("taxonomy is a required parameter for the classify.seqs command."); 
				mothurOutEndLine();
				abort = true; 
			}
			else if (taxonomyFileName == "not open") { abort = true; }	

			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "ksize", false);		if (temp == "not found"){	temp = "8";				}
			convert(temp, kmerSize); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
			search = validParameter.validFile(parameters, "search", false);		if (search == "not found"){	search = "kmer";		}
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found"){	method = "bayesian";	}
			
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			convert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, misMatch);  
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			convert(temp, gapOpen);  
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, gapExtend); 
			
			temp = validParameter.validFile(parameters, "numwanted", false);	if (temp == "not found"){	temp = "10";			}
			convert(temp, numWanted);
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found"){	temp = "0";				}
			convert(temp, cutoff);
			
			temp = validParameter.validFile(parameters, "probs", false);		if (temp == "not found"){	temp = "true";			}
			probs = isTrue(temp);

			
			if ((method == "bayesian") && (search != "kmer"))  { 
				mothurOut("The bayesian method requires the kmer search." + search + "will be disregarded." ); mothurOutEndLine();
				search = "kmer";
			}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "ClassifySeqsCommand", "ClassifySeqsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

ClassifySeqsCommand::~ClassifySeqsCommand(){	

	if (abort == false) {
		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
	}
}

//**********************************************************************************************************************

void ClassifySeqsCommand::help(){
	try {
		mothurOut("The classify.seqs command reads a fasta file containing sequences and creates a .taxonomy file and a .tax.summary file.\n");
		mothurOut("The classify.seqs command parameters are template, fasta, search, ksize, method, taxonomy, processors, match, mismatch, gapopen, gapextend, numwanted and probs.\n");
		mothurOut("The template, fasta and taxonomy parameters are required.\n");
		mothurOut("The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer.\n");
		mothurOut("The method parameter allows you to specify classification method to use.  Your options are: bayesian and knn. The default is bayesian.\n");
		mothurOut("The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 8.\n");
		mothurOut("The processors parameter allows you to specify the number of processors to use. The default is 1.\n");
		mothurOut("The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n");
		mothurOut("The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n");
		mothurOut("The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -1.0.\n");
		mothurOut("The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -2.0.\n");
		mothurOut("The numwanted parameter allows you to specify the number of sequence matches you want with the knn method.  The default is 10.\n");
		mothurOut("The cutoff parameter allows you to specify a bootstrap confidence threshold for your taxonomy.  The default is 0.\n");
		mothurOut("The probs parameter shut off the bootstrapping results for the bayesian method. The default is true, meaning you want the bootstrapping to be run.\n");
		mothurOut("The classify.seqs command should be in the following format: \n");
		mothurOut("classify.seqs(template=yourTemplateFile, fasta=yourFastaFile, method=yourClassificationMethod, search=yourSearchmethod, ksize=yourKmerSize, taxonomy=yourTaxonomyFile, processors=yourProcessors) \n");
		mothurOut("Example classify.seqs(fasta=amazon.fasta, template=core.filtered, method=knn, search=gotoh, ksize=8, processors=2)\n");
		mothurOut("The .taxonomy file consists of 2 columns: 1 = your sequence name, 2 = the taxonomy for your sequence. \n");
		mothurOut("The .tax.summary is a summary of the different taxonomies represented in your fasta file. \n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "ClassifySeqsCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

int ClassifySeqsCommand::execute(){
	try {
		if (abort == true) {	return 0;	}
		
		if(method == "bayesian")			{	classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, probs);		}
		else if(method == "knn")			{	classify = new Knn(taxonomyFileName, templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch, numWanted);				}
		else {
			mothurOut(search + " is not a valid method option. I will run the command using bayesian.");
			mothurOutEndLine();
			classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, probs);	
		}

		int numFastaSeqs = 0;
		
		string newTaxonomyFile = getRootName(fastaFileName) + "taxonomy";
		string tempTaxonomyFile = getRootName(fastaFileName) + "taxonomy.temp";
		string taxSummary = getRootName(fastaFileName) + "tax.summary";
		
		int start = time(NULL);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors == 1){
			ifstream inFASTA;
			openInputFile(fastaFileName, inFASTA);
			numFastaSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
			inFASTA.close();
			
			lines.push_back(new linePair(0, numFastaSeqs));
		
			driver(lines[0], newTaxonomyFile, tempTaxonomyFile);
		}
		else{
			vector<int> positions;
			processIDS.resize(0);
			
			ifstream inFASTA;
			openInputFile(fastaFileName, inFASTA);
			
			string input;
			while(!inFASTA.eof()){
				input = getline(inFASTA);
				if (input.length() != 0) {
					if(input[0] == '>'){	int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
				}
			}
			inFASTA.close();
			
			numFastaSeqs = positions.size();
			
			int numSeqsPerProcessor = numFastaSeqs / processors;
			
			for (int i = 0; i < processors; i++) {
				int startPos = positions[ i * numSeqsPerProcessor ];
				if(i == processors - 1){
					numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;
				}
				lines.push_back(new linePair(startPos, numSeqsPerProcessor));
			}
			createProcesses(newTaxonomyFile, tempTaxonomyFile); 
			
			rename((newTaxonomyFile + toString(processIDS[0]) + ".temp").c_str(), newTaxonomyFile.c_str());
			rename((tempTaxonomyFile + toString(processIDS[0]) + ".temp").c_str(), tempTaxonomyFile.c_str());
			
			for(int i=1;i<processors;i++){
				appendTaxFiles((newTaxonomyFile + toString(processIDS[i]) + ".temp"), newTaxonomyFile);
				appendTaxFiles((tempTaxonomyFile + toString(processIDS[i]) + ".temp"), tempTaxonomyFile);
				remove((newTaxonomyFile + toString(processIDS[i]) + ".temp").c_str());
				remove((tempTaxonomyFile + toString(processIDS[i]) + ".temp").c_str());
			}
			
		}
#else
		ifstream inFASTA;
		openInputFile(fastaFileName, inFASTA);
		numFastaSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
		inFASTA.close();
		
		lines.push_back(new linePair(0, numFastaSeqs));
		
		driver(lines[0], newTaxonomyFile, tempTaxonomyFile);
#endif	
		delete classify;
		
		//make taxonomy tree from new taxonomy file 
		ifstream inTaxonomy;
		openInputFile(tempTaxonomyFile, inTaxonomy);
		
		string accession, taxaList;
		PhyloTree taxaBrowser;
		
		//read in users taxonomy file and add sequences to tree
		while(!inTaxonomy.eof()){
			inTaxonomy >> accession >> taxaList;
			
			taxaBrowser.addSeqToTree(accession, taxaList);
			
			gobble(inTaxonomy);
		}
		inTaxonomy.close();
		remove(tempTaxonomyFile.c_str());
		
		taxaBrowser.assignHeirarchyIDs(0);
		
		ofstream outTaxTree;
		openOutputFile(taxSummary, outTaxTree);
		
		taxaBrowser.print(outTaxTree);
		
		mothurOutEndLine();
		mothurOut("It took " + toString(time(NULL) - start) + " secs to classify " + toString(numFastaSeqs) + " sequences.");
		mothurOutEndLine();
		mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "ClassifySeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/

void ClassifySeqsCommand::createProcesses(string taxFileName, string tempTaxFile) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		//		processIDS.resize(0);
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driver(lines[process], taxFileName + toString(getpid()) + ".temp", tempTaxFile + toString(getpid()) + ".temp");
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
#endif		
	}
	catch(exception& e) {
		errorOut(e, "ClassifySeqsCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/

void ClassifySeqsCommand::appendTaxFiles(string temp, string filename) {
	try{
		
		ofstream output;
		ifstream input;
		openOutputFileAppend(filename, output);
		openInputFile(temp, input);
		
		while(char c = input.get()){
			if(input.eof())		{	break;			}
			else				{	output << c;	}
		}
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		errorOut(e, "ClassifySeqsCommand", "appendTaxFiles");
		exit(1);
	}
}

//**********************************************************************************************************************

int ClassifySeqsCommand::driver(linePair* line, string taxFName, string tempTFName){
	try {
		ofstream outTax;
		openOutputFile(taxFName, outTax);
		
		ofstream outTaxSimple;
		openOutputFile(tempTFName, outTaxSimple);
	
		ifstream inFASTA;
		openInputFile(fastaFileName, inFASTA);

		inFASTA.seekg(line->start);
		
		string taxonomy;

		for(int i=0;i<line->numSeqs;i++){
			
			Sequence* candidateSeq = new Sequence(inFASTA);
			
			if (candidateSeq->getName() != "") {
				taxonomy = classify->getTaxonomy(candidateSeq);
				
				if (taxonomy != "bad seq") {
					outTax << candidateSeq->getName() << '\t' << taxonomy << endl;
					outTaxSimple << candidateSeq->getName() << '\t' << classify->getSimpleTax() << endl;
				}
			}				
			delete candidateSeq;
			
			if((i+1) % 100 == 0){
				mothurOut("Classifying sequence " + toString(i+1)); mothurOutEndLine();
			}
		}
		
		inFASTA.close();
		outTax.close();
		outTaxSimple.close();
		
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "ClassifySeqsCommand", "driver");
		exit(1);
	}
}

/**************************************************************************************************/
