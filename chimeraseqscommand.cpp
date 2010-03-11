/*
 *  chimeraseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimeraseqscommand.h"
#include "bellerophon.h"
#include "pintail.h"
#include "ccode.h"
#include "chimeracheckrdp.h"
#include "chimeraslayer.h"


//***************************************************************************************************************

ChimeraSeqsCommand::ChimeraSeqsCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "filter", "correction", "processors", "method", "window", "increment", "template", "conservation", "quantile", "mask", 
			"numwanted", "ksize", "svg", "name", "match","mismatch", "divergence", "minsim","mincov","minbs", "minsnp","parents", "iters","outputdir","inputdir", "search","realign" };
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("template");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["template"] = inputDir + it->second;		}
				}
				
				it = parameters.find("conservation");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["conservation"] = inputDir + it->second;		}
				}
				
				it = parameters.find("quantile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["quantile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the chimera.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			templatefile = validParameter.validFile(parameters, "template", true);
			if (templatefile == "not open") { abort = true; }
			else if (templatefile == "not found") { templatefile = "";  }	
			
			consfile = validParameter.validFile(parameters, "conservation", true);
			if (consfile == "not open") { abort = true; }
			else if (consfile == "not found") { consfile = "";  }	
			
			quanfile = validParameter.validFile(parameters, "quantile", true);
			if (quanfile == "not open") { abort = true; }
			else if (quanfile == "not found") { quanfile = "";  }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = "";  }

			maskfile = validParameter.validFile(parameters, "mask", false);
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  { 
				if (inputDir != "") {
					string path = hasPath(maskfile);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	maskfile = inputDir + maskfile;		}
				}

				ifstream in;
				int	ableToOpen = openInputFile(maskfile, in);
				if (ableToOpen == 1) { abort = true; }
				in.close();
			}
			
			method = validParameter.validFile(parameters, "method", false);			if (method == "not found") { method = "pintail"; }
			
			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "F"; }
			filter = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "correction", false);		if (temp == "not found") { temp = "T"; }
			correction = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "7"; }
			convert(temp, ksize);
			
			temp = validParameter.validFile(parameters, "svg", false);				if (temp == "not found") { temp = "F"; }
			svg = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "window", false);	
			if ((temp == "not found") && (method == "chimeraslayer")) { temp = "50"; }			
			else if (temp == "not found") { temp = "0"; }
			convert(temp, window);
			
			temp = validParameter.validFile(parameters, "match", false);			if (temp == "not found") { temp = "5"; }
			convert(temp, match);
			
			temp = validParameter.validFile(parameters, "mismatch", false);			if (temp == "not found") { temp = "-4"; }
			convert(temp, mismatch);
			
			temp = validParameter.validFile(parameters, "divergence", false);		if (temp == "not found") { temp = "1.007"; }
			convert(temp, divR);
			
			temp = validParameter.validFile(parameters, "minsim", false);			if (temp == "not found") { temp = "90"; }
			convert(temp, minSimilarity);
			
			temp = validParameter.validFile(parameters, "mincov", false);			if (temp == "not found") { temp = "70"; }
			convert(temp, minCoverage);
			
			temp = validParameter.validFile(parameters, "minbs", false);			if (temp == "not found") { temp = "90"; }
			convert(temp, minBS);
			
			temp = validParameter.validFile(parameters, "minsnp", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, minSNP);

			temp = validParameter.validFile(parameters, "parents", false);			if (temp == "not found") { temp = "3"; }
			convert(temp, parents); 
			
			temp = validParameter.validFile(parameters, "realign", false);			if (temp == "not found") { temp = "f"; }
			realign = isTrue(temp); 
			
			search = validParameter.validFile(parameters, "search", false);			if (search == "not found") { search = "distance"; }
			
			temp = validParameter.validFile(parameters, "iters", false);	
			if ((temp == "not found") && (method == "chimeraslayer")) { temp = "100"; }		
			else if (temp == "not found") { temp = "1000"; }
			convert(temp, iters); 
			 
			temp = validParameter.validFile(parameters, "increment", false);		
			if ((temp == "not found") && (method == "chimeracheck")) { temp = "10"; }
			else if ((temp == "not found") && (method == "chimeraslayer")) { temp = "5"; }
			else if (temp == "not found") { temp = "25"; }
			convert(temp, increment);
			
			temp = validParameter.validFile(parameters, "numwanted", false);
			if ((temp == "not found") && (method == "chimeraslayer")) { temp = "15"; }		
			else if (temp == "not found") { temp = "20"; }
			convert(temp, numwanted);

			if ((search != "distance") && (search != "blast") && (search != "kmer")) { m->mothurOut(search + " is not a valid search."); m->mothurOutEndLine(); abort = true;  }
			
			if (((method != "bellerophon")) && (templatefile == "")) { m->mothurOut("You must provide a template file with the pintail, ccode, chimeraslayer or chimeracheck methods."); m->mothurOutEndLine(); abort = true;  }
			

		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "ChimeraSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraSeqsCommand::help(){
	try {
	
		//"fasta", "filter", "correction", "processors", "method", "window", "increment", "template", "conservation", "quantile", "mask", "numwanted", "ksize", "svg", "name"
		//m->mothurOut("chimera.seqs ASSUMES that your sequences are ALIGNED and if using a template that the template file sequences are the same length as the fasta file sequences.\n\n");
		m->mothurOut("The chimera.seqs command reads a fastafile and creates list of potentially chimeric sequences.\n");
		m->mothurOut("The chimera.seqs command parameters are fasta, filter, correction, processors, mask, method, window, increment, template, conservation, quantile, numwanted, ksize, svg, name, iters, search, realign.\n");
		m->mothurOut("The fasta parameter is always required and template is required if using pintail, ccode or chimeracheck.\n");
		m->mothurOut("The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter. \n");
		m->mothurOut("The correction parameter allows you to put more emphasis on the distance between highly similar sequences and less emphasis on the differences between remote homologs.\n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		m->mothurOut("The method parameter allows you to specify the method for finding chimeric sequences.  The default is pintail. Options include bellerophon, ccode and chimeracheck \n");
		m->mothurOut("The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences. \n");
		m->mothurOut("The window parameter allows you to specify the window size for searching for chimeras. \n");
		m->mothurOut("The increment parameter allows you to specify how far you move each window while finding chimeric sequences.\n");
		m->mothurOut("The template parameter allows you to enter a template file containing known non-chimeric sequences. \n");
		m->mothurOut("The conservation parameter allows you to enter a frequency file containing the highest bases frequency at each place in the alignment.\n");
		m->mothurOut("The quantile parameter allows you to enter a file containing quantiles for a template files sequences.\n");
		m->mothurOut("The numwanted parameter allows you to specify how many sequences you would each query sequence compared with.\n");
		m->mothurOut("The ksize parameter allows you to input kmersize. \n");
		m->mothurOut("The svg parameter allows you to specify whether or not you would like a svg file outputted for each query sequence.\n");
		m->mothurOut("The name parameter allows you to enter a file containing names of sequences you would like .svg files for.\n");
		m->mothurOut("The iters parameter allows you to specify the number of bootstrap iters to do with the chimeraslayer method.\n");
		m->mothurOut("The minsim parameter allows you .... \n");
		m->mothurOut("The mincov parameter allows you to specify minimum coverage by closest matches found in template. Default is 70, meaning 70%. \n");
		m->mothurOut("The minbs parameter allows you to specify minimum bootstrap support for calling a sequence chimeric. Default is 90, meaning 90%. \n");
		m->mothurOut("The minsnp parameter allows you to specify percent of SNPs to sample on each side of breakpoint for computing bootstrap support (default: 10) \n");
		m->mothurOut("The search parameter allows you to specify search method for finding the closest parent. Choices are distance, blast, and kmer, default distance.  -used only by chimeraslayer. \n");
		m->mothurOut("The realign parameter allows you to realign the query to the potential paretns. Choices are true or false, default false.  -used only by chimeraslayer. \n");
		m->mothurOut("NOT ALL PARAMETERS ARE USED BY ALL METHODS. Please look below for method specifics.\n\n");
		m->mothurOut("Details for each method: \n"); 
		m->mothurOut("\tpintail: \n"); 
		m->mothurOut("\t\tparameters: fasta=required, template=required, filter=F, mask=no mask, processors=1, window=300, increment=25, conservation=not required, but will improve speed, quantile=not required, but will greatly improve speed. \n"); 
		m->mothurOut("\t\tIf you have run chimera.seqs using pintail a .quan and .freq file will be created for your template, if you have not provided them for use in future command executions.\n");
		m->mothurOut("\tbellerophon: \n"); 
		m->mothurOut("\t\tparameters: fasta=required, filter=F, processors=1, window=1/4 length of seq, increment=25, correction=T. \n"); 
		m->mothurOut("\tccode: \n"); 
		m->mothurOut("\t\tparameters: fasta=required, template=required, filter=F, mask=no mask, processors=1, window=10% of length, numwanted=20\n"); 
		m->mothurOut("\tchimeracheck: \n"); 
		m->mothurOut("\t\tparameters: fasta=required, template=required, processors=1, increment=10, ksize=7, svg=F, name=none\n\n"); 
		m->mothurOut("\tchimeraslayer: \n"); 
		m->mothurOut("\t\tparameters: fasta=required, template=required, processors=1, increment=10, mask=no mask, numwanted=10, match=5, mismatch=-4, divergence=1.0, minsim=90, parents=5, iters=1000, window=100. \n\n"); 
		m->mothurOut("The chimera.seqs command should be in the following format: \n");
		m->mothurOut("chimera.seqs(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors, method=bellerophon) \n");
		m->mothurOut("Example: chimera.seqs(fasta=AD.align, filter=True, correction=true, method=bellerophon, window=200) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraSeqsCommand::~ChimeraSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		int start = time(NULL);	
		
		if (method == "bellerophon")			{		chimera = new Bellerophon(fastafile, outputDir);			}
		else if (method == "pintail")			{		chimera = new Pintail(fastafile, outputDir);				}
		else if (method == "ccode")				{		chimera = new Ccode(fastafile, outputDir);					}
		else if (method == "chimeracheck")		{		chimera = new ChimeraCheckRDP(fastafile, outputDir);		}
		else if (method == "chimeraslayer")		{		chimera = new ChimeraSlayer(search, realign, fastafile);	}
		else { m->mothurOut("Not a valid method."); m->mothurOutEndLine(); return 0;		}
		
		//set user options
		if (maskfile == "default") { m->mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); m->mothurOutEndLine();  }
		
		chimera->setCons(consfile);	
		chimera->setQuantiles(quanfile);				
		chimera->setMask(maskfile);
		chimera->setFilter(filter);
		chimera->setCorrection(correction);
		chimera->setProcessors(processors);
		chimera->setWindow(window);
		chimera->setIncrement(increment);
		chimera->setNumWanted(numwanted);
		chimera->setKmerSize(ksize);
		chimera->setSVG(svg);
		chimera->setName(namefile);
		chimera->setMatch(match);
		chimera->setMisMatch(mismatch);
		chimera->setDivR(divR);
		chimera->setParents(parents);
		chimera->setMinSim(minSimilarity);
		chimera->setMinCoverage(minCoverage);
		chimera->setMinBS(minBS);
		chimera->setMinSNP(minSNP);
		chimera->setIters(iters);
		

		string outputFileName = outputDir + getRootName(getSimpleName(fastafile)) + method + maskfile + ".chimeras";
		string accnosFileName = outputDir + getRootName(getSimpleName(fastafile)) + method + maskfile + ".accnos";
		bool hasAccnos = true;
		
		if (method == "bellerophon") {//run bellerophon separately since you need to read entire fastafile to run it
			chimera->getChimeras();
			
			if (m->control_pressed) { delete chimera;	return 0;	}
			
			ofstream out;
			openOutputFile(outputFileName, out);
			
			ofstream out2;
			openOutputFile(accnosFileName, out2);
			
			chimera->print(out, out2);
			out.close();
			out2.close(); 
			
			if (m->control_pressed) { remove(accnosFileName.c_str()); remove(outputFileName.c_str()); delete chimera;	return 0;	}
			
			//delete accnos file if its blank 
			if (isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
			
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			m->mothurOut(outputFileName); m->mothurOutEndLine();	
			if (hasAccnos) {  m->mothurOut(accnosFileName); m->mothurOutEndLine();  }
			m->mothurOutEndLine();
			
			delete chimera;
			return 0;
		}
		
		//reads template
		chimera->setTemplateFile(templatefile);
		
		if (m->control_pressed) { delete chimera;	return 0;	}
		
		if  (method != "chimeracheck") {   
			if (chimera->getUnaligned()) { 
				m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); 
				delete chimera;
				return 0; 
			}
		}
		
		//some methods need to do prep work before processing the chimeras
		chimera->doPrep(); 
		
		if (m->control_pressed) { delete chimera;	return 0;	}
		
		templateSeqsLength = chimera->getLength();
		
		ofstream outHeader;
		string tempHeader = outputDir + getRootName(getSimpleName(fastafile)) + method + maskfile + ".chimeras.tempHeader";
		openOutputFile(tempHeader, outHeader);
		
		chimera->printHeader(outHeader);
		outHeader.close();
		
		
		//break up file
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				ifstream inFASTA;
				openInputFile(fastafile, inFASTA);
				numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
				inFASTA.close();
				
				lines.push_back(new linePair(0, numSeqs));
				
				driver(lines[0], outputFileName, fastafile, accnosFileName);
				
				if (m->control_pressed) { 
					remove(outputFileName.c_str()); 
					remove(tempHeader.c_str()); 
					remove(accnosFileName.c_str());
					for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
					delete chimera;
					return 0;
				}
				
				//delete accnos file if its blank 
				if (isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
								
			}else{
				vector<int> positions;
				processIDS.resize(0);
				
				ifstream inFASTA;
				openInputFile(fastafile, inFASTA);
				
				string input;
				while(!inFASTA.eof()){
					input = getline(inFASTA);
					if (input.length() != 0) {
						if(input[0] == '>'){	long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
					}
				}
				inFASTA.close();
				
				numSeqs = positions.size();
				
				int numSeqsPerProcessor = numSeqs / processors;
				
				for (int i = 0; i < processors; i++) {
					long int startPos = positions[ i * numSeqsPerProcessor ];
					if(i == processors - 1){
						numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor;
					}
					lines.push_back(new linePair(startPos, numSeqsPerProcessor));
				}
				
				
				createProcesses(outputFileName, fastafile, accnosFileName); 
			
				rename((outputFileName + toString(processIDS[0]) + ".temp").c_str(), outputFileName.c_str());
					
				//append output files
				for(int i=1;i<processors;i++){
					appendOutputFiles((outputFileName + toString(processIDS[i]) + ".temp"), outputFileName);
					remove((outputFileName + toString(processIDS[i]) + ".temp").c_str());
				}
				
				vector<string> nonBlankAccnosFiles;
				//delete blank accnos files generated with multiple processes
				for(int i=0;i<processors;i++){  
					if (!(isBlank(accnosFileName + toString(processIDS[i]) + ".temp"))) {
						nonBlankAccnosFiles.push_back(accnosFileName + toString(processIDS[i]) + ".temp");
					}else { remove((accnosFileName + toString(processIDS[i]) + ".temp").c_str());  }
				}
				
				//append accnos files
				if (nonBlankAccnosFiles.size() != 0) { 
					rename(nonBlankAccnosFiles[0].c_str(), accnosFileName.c_str());
					
					for (int h=1; h < nonBlankAccnosFiles.size(); h++) {
						appendOutputFiles(nonBlankAccnosFiles[h], accnosFileName);
						remove(nonBlankAccnosFiles[h].c_str());
					}
				}else{ hasAccnos = false;  }
				
				if (m->control_pressed) { 
					remove(outputFileName.c_str()); 
					remove(accnosFileName.c_str());
					for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
					delete chimera;
					return 0;
				}

			}

		#else
			ifstream inFASTA;
			openInputFile(candidateFileNames[s], inFASTA);
			numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
			inFASTA.close();
			lines.push_back(new linePair(0, numSeqs));
			
			driver(lines[0], outputFileName, fastafile, accnosFileName);
			
			if (m->control_pressed) { 
					remove(outputFileName.c_str()); 
					remove(tempHeader.c_str()); 
					remove(accnosFileName.c_str());
					for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
					delete chimera;
					return 0;
			}
			
			//delete accnos file if its blank 
			if (isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
		#endif
		
		//m->mothurOut("Output File Names: ");
		//if ((filter) && (method == "bellerophon")) { m->mothurOut(
		//if (outputDir == "") { fastafile = getRootName(fastafile) + "filter.fasta"; }
		//	else				 { fastafile = outputDir + getRootName(getSimpleName(fastafile)) + "filter.fasta"; }
	
		appendOutputFiles(tempHeader, outputFileName);
	
		remove(outputFileName.c_str());
		rename(tempHeader.c_str(), outputFileName.c_str());
	
		delete chimera;
		
		if (method == "chimeracheck") { remove(accnosFileName.c_str());  m->mothurOutEndLine(); m->mothurOut("This method does not determine if a sequence is chimeric, but allows you to make that determination based on the IS values."); m->mothurOutEndLine();  }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		if (hasAccnos) {  m->mothurOut(accnosFileName); m->mothurOutEndLine();  }
		m->mothurOutEndLine();

		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
		
		m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "execute");
		exit(1);
	}
}//**********************************************************************************************************************

int ChimeraSeqsCommand::driver(linePair* line, string outputFName, string filename, string accnos){
	try {
		ofstream out;
		openOutputFile(outputFName, out);
		
		ofstream out2;
		openOutputFile(accnos, out2);
		
		ifstream inFASTA;
		openInputFile(filename, inFASTA);

		inFASTA.seekg(line->start);
		
		for(int i=0;i<line->numSeqs;i++){
		
			if (m->control_pressed) {	return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if ((candidateSeq->getAligned().length() != templateSeqsLength) && (method != "chimeracheck")) {  //chimeracheck does not require seqs to be aligned
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (m->control_pressed) {	delete candidateSeq; return 1;	}
		
					//print results
					chimera->print(out, out2);
				}
			}
			delete candidateSeq;
			
			//report progress
			if((i+1) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(i+1)); m->mothurOutEndLine();		}
		}
		//report progress
		if((line->numSeqs) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(line->numSeqs)); m->mothurOutEndLine();		}
		
		out.close();
		out2.close();
		inFASTA.close();
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "driver");
		exit(1);
	}
}

/**************************************************************************************************/

int ChimeraSeqsCommand::createProcesses(string outputFileName, string filename, string accnos) {
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
				driver(lines[process], outputFileName + toString(getpid()) + ".temp", filename, accnos + toString(getpid()) + ".temp");
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		return 0;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/

void ChimeraSeqsCommand::appendOutputFiles(string temp, string filename) {
	try{
		
		ofstream output;
		ifstream input;
		
		openOutputFileAppend(temp, output);
		openInputFile(filename, input, "noerror");
		
		while(char c = input.get()){
			if(input.eof())		{	break;			}
			else				{	output << c;	}
		}
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSeqsCommand", "appendOuputFiles");
		exit(1);
	}
}
//**********************************************************************************************************************


