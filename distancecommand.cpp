/*
 *  distancecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "distancecommand.h"
#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"

//**********************************************************************************************************************

DistanceCommand::DistanceCommand(string option) {
	try {
		abort = false;
		Estimators.clear();
				
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "output", "calc", "countends", "cutoff", "processors", "outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it2;
		
			//check to make sure all parameters are valid for command
			for (it2 = parameters.begin(); it2 != parameters.end(); it2++) { 
				if (validParameter.isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it2 = parameters.find("fasta");
				//user has given a template file
				if(it2 != parameters.end()){ 
					path = hasPath(it2->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it2->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { m->mothurOut("fasta is a required parameter for the dist.seqs command."); m->mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
			else{
				ifstream inFASTA;
				openInputFile(fastafile, inFASTA);
				alignDB = SequenceDB(inFASTA); 
				inFASTA.close();
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "onegap";  }
			else { 
				 if (calc == "default")  {  calc = "onegap";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "countends", false);	if(temp == "not found"){	temp = "T";	}
			convert(temp, countends); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if(temp == "not found"){	temp = "1.0"; }
			convert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if(temp == "not found"){	temp = "1"; }
			convert(temp, processors); 
			
			output = validParameter.validFile(parameters, "output", false);		if(output == "not found"){	output = "column"; }
			
			if ((output != "column") && (output != "lt") && (output != "square")) { m->mothurOut(output + " is not a valid output form. Options are column, lt and square. I will use column."); m->mothurOutEndLine(); output = "column"; }
			
			ValidCalculators validCalculator;
			
			if (isTrue(countends) == true) {
				for (int i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
						if (Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
						else if (Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
						else if (Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
					}
				}
			}else {
				for (int i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
						if (Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
						else if (Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
						else if (Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
					}
				}
			}

		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "DistanceCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

DistanceCommand::~DistanceCommand(){
	
	for(int i=0;i<lines.size();i++){
		delete lines[i];
	}
	
}
	
//**********************************************************************************************************************

void DistanceCommand::help(){
	try {
		m->mothurOut("The dist.seqs command reads a file containing sequences and creates a distance file.\n");
		m->mothurOut("The dist.seqs command parameters are fasta, calc, countends, output, cutoff and processors.  \n");
		m->mothurOut("The fasta parameter is required.\n");
		m->mothurOut("The calc parameter allows you to specify the method of calculating the distances.  Your options are: nogaps, onegap or eachgap. The default is onegap.\n");
		m->mothurOut("The countends parameter allows you to specify whether to include terminal gaps in distance.  Your options are: T or F. The default is T.\n");
		m->mothurOut("The cutoff parameter allows you to specify maximum distance to keep. The default is 1.0.\n");
		m->mothurOut("The output parameter allows you to specify format of your distance matrix. Options are column, lt, and square. The default is column.\n");
		m->mothurOut("The processors parameter allows you to specify number of processors to use.  The default is 1.\n");
		m->mothurOut("The dist.seqs command should be in the following format: \n");
		m->mothurOut("dist.seqs(fasta=yourFastaFile, calc=yourCalc, countends=yourEnds, cutoff= yourCutOff, processors=yourProcessors) \n");
		m->mothurOut("Example dist.seqs(fasta=amazon.fasta, calc=eachgap, countends=F, cutoff= 2.0, processors=3).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. calc), '=' and parameters (i.e.yourCalc).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************

int DistanceCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		int numSeqs = alignDB.getNumSeqs();
		cutoff += 0.005;
		
		string outputFile;
				
		if (output == "lt") { //does the user want lower triangle phylip formatted file 
			outputFile = outputDir + getRootName(getSimpleName(fastafile)) + "phylip.dist";
			remove(outputFile.c_str());
			
			//output numSeqs to phylip formatted dist file
		}else if (output == "column") { //user wants column format
			outputFile = outputDir + getRootName(getSimpleName(fastafile)) + "dist";
			remove(outputFile.c_str());
		}else { //assume square
			outputFile = outputDir + getRootName(getSimpleName(fastafile)) + "square.dist";
			remove(outputFile.c_str());
		}
		

#ifdef USE_MPI
		
		int pid, start, end; 
		int tag = 2001;
				
		MPI_Status status; 
		MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		//each process gets where it should start and stop in the file
		start = int (sqrt(float(pid)/float(processors)) * numSeqs);
		end = int (sqrt(float(pid+1)/float(processors)) * numSeqs);
		
		if (pid == 0) { //you are the root process 
			//do your part
			string outputMyPart;
			driverMPI(start, end, outputMyPart, cutoff);
			
			ofstream out;
			openOutputFile(outputFile, out);
			
			out << outputMyPart;
			
			//get the childrens parts
			for(int i = 1; i < processors; i++) { 
				int length;
				MPI_Recv(&length, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status); 
				
				char buf[length];
					
				MPI_Recv(buf, length, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status); 
				
				outputMyPart = buf;
				out << outputMyPart;
			}
			
			out.close();
			
		}else { //you are a child process
			//do your part
			string outputMyPart;
			driverMPI(start, end, outputMyPart, cutoff);
		
			//send results to parent
			int length = outputMyPart.length();
			char buf[length];
			strcpy(buf, outputMyPart.c_str()); 
			
			MPI_Send( &length, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(buf, length, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
		}
		

#else		
				
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		//if you don't need to fork anything
		if(processors == 1){
			driver(0, numSeqs, outputFile, cutoff);
		}else{ //you have multiple processors
			
			for (int i = 0; i < processors; i++) {
				lines.push_back(new linePair());
				lines[i]->start = int (sqrt(float(i)/float(processors)) * numSeqs);
				lines[i]->end = int (sqrt(float(i+1)/float(processors)) * numSeqs);
			}

			createProcesses(outputFile); 
		
			map<int, int>::iterator it = processIDS.begin();
			rename((outputFile + toString(it->second) + ".temp").c_str(), outputFile.c_str());
			it++;
			
			//append and remove temp files
			for (; it != processIDS.end(); it++) {
				appendFiles((outputFile + toString(it->second) + ".temp"), outputFile);
				remove((outputFile + toString(it->second) + ".temp").c_str());
			}
		}
	#else
		ifstream inFASTA;
		driver(0, numSeqs, outputFile, cutoff);
	#endif
	
#endif
		if (m->control_pressed) { delete distCalculator; remove(outputFile.c_str()); return 0; }
		
		#ifdef USE_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		if (output == "square") {  convertMatrix(outputFile); }
		
		#ifdef USE_MPI
			}
		#endif
		
		if (m->control_pressed) { delete distCalculator; remove(outputFile.c_str()); return 0; }
		
		delete distCalculator;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFile); m->mothurOutEndLine();
		m->mothurOutEndLine();

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
void DistanceCommand::createProcesses(string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS[lines[process]->end] = pid;  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driver(lines[process]->start, lines[process]->end, filename + toString(getpid()) + ".temp", cutoff);
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
	
		//force parent to wait until all the processes are done
		for (map<int, int>::iterator it = processIDS.begin(); it != processIDS.end(); it++) { 
			int temp = it->second;
			wait(&temp);
		}
#endif
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driver(int startLine, int endLine, string dFileName, float cutoff){
	try {

		int startTime = time(NULL);
		
		//column file
		ofstream outFile(dFileName.c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		
		if((output == "lt") && startLine == 0){	outFile << alignDB.getNumSeqs() << endl;	}
		
		for(int i=startLine;i<endLine;i++){
			if(output == "lt")	{	
				string name = alignDB.get(i).getName();
				if (name.length() < 10) { //pad with spaces to make compatible
					while (name.length() < 10) {  name += " ";  }
				}
				outFile << name << '\t';	
			}
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) { outFile.close(); return 0;  }
				
				distCalculator->calcDist(alignDB.get(i), alignDB.get(j));
				double dist = distCalculator->getDist();
				
				if(dist <= cutoff){
					if (output == "column") { outFile << alignDB.get(i).getName() << ' ' << alignDB.get(j).getName() << ' ' << dist << endl; }
				}
				if (output == "lt") {  outFile << dist << '\t'; }
				
				if (output == "square") { //make a square column you can convert to square phylip
					outFile << alignDB.get(i).getName() << '\t' << alignDB.get(j).getName() << '\t' << dist << endl;
					outFile << alignDB.get(j).getName() << '\t' << alignDB.get(i).getName() << '\t' << dist << endl;
				}

			}
			
			if (output == "lt") { outFile << endl; }
			
			if(i % 100 == 0){
				m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
			}
			
		}
		m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		
		outFile.close();
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driverMPI(int startLine, int endLine, string& outputString, float cutoff){
	try {

		int startTime = time(NULL);
		
		outputString = "";
				
		if((output == "lt") && startLine == 0){	outputString += (toString(alignDB.getNumSeqs()) + '\n');	}
		
		for(int i=startLine;i<endLine;i++){
	
			if(output == "lt")	{	
				string name = alignDB.get(i).getName();
				if (name.length() < 10) { //pad with spaces to make compatible
					while (name.length() < 10) {  name += " ";  }
				}
				outputString += (name + '\t');	
			}
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) {  return 0;  }
				
				distCalculator->calcDist(alignDB.get(i), alignDB.get(j));
				double dist = distCalculator->getDist();
				
				if(dist <= cutoff){
					if (output == "column") { outputString += (alignDB.get(i).getName() + ' ' + alignDB.get(j).getName() + ' ' + toString(dist) + '\n'); }
				}
				if (output == "lt") {   outputString += (toString(dist) + '\t'); }
				
				if (output == "square") { //make a square column you can convert to square phylip
					outputString += (alignDB.get(i).getName() + ' ' + alignDB.get(j).getName() + ' ' + toString(dist) + '\n');
					outputString += (alignDB.get(j).getName() + ' ' + alignDB.get(i).getName() + ' ' + toString(dist) + '\n');
				}

			}
			
			if (output == "lt") { outputString += '\n'; }
			
			if(i % 100 == 0){
				m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
			}
			
		}
		m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "driver");
		exit(1);
	}
}

/**************************************************************************************************/
int DistanceCommand::convertMatrix(string outputFile) {
	try{

		//sort file by first column so the distances for each row are together
		string outfile = getRootName(outputFile) + "sorted.dist.temp";
		
		//use the unix sort 
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			string command = "sort -n " + outputFile + " -o " + outfile;
			system(command.c_str());
		#else //sort using windows sort
			string command = "sort " + outputFile + " /O " + outfile;
			system(command.c_str());
		#endif
		

		//output to new file distance for each row and save positions in file where new row begins
		ifstream in;
		openInputFile(outfile, in);
		
		ofstream out;
		openOutputFile(outputFile, out);
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);

		out << alignDB.getNumSeqs() << endl;
		
		//get first currentRow
		string first, currentRow, second;
		float dist;
		map<string, float> rowDists; //take advantage of the fact that maps are already sorted by key 
		map<string, float>::iterator it;
		
		in >> first;
		currentRow = first;
		
		rowDists[first] = 0.00; //distance to yourself is 0.0
		
		in.seekg(0);
		//openInputFile(outfile, in);
		
		while(!in.eof()) {
			if (m->control_pressed) { in.close(); remove(outfile.c_str()); out.close(); return 0; }
			
			in >> first >> second >> dist; gobble(in);
				
			if (first != currentRow) {
				//print out last row
				out << currentRow << '\t'; //print name

				//print dists
				for (it = rowDists.begin(); it != rowDists.end(); it++) {
					out << it->second << '\t';
				}
				out << endl;
				
				//start new row
				currentRow = first;
				rowDists.clear();
				rowDists[first] = 0.00;
				rowDists[second] = dist;
			}else{
				rowDists[second] = dist;
			}
		}
		//print out last row
		out << currentRow << '\t'; //print name
				
		//print dists
		for (it = rowDists.begin(); it != rowDists.end(); it++) {
			out << it->second << '\t';
		}
		out << endl;
		
		in.close();
		out.close();
		
		remove(outfile.c_str());
		
		return 1;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "convertMatrix");
		exit(1);
	}
}
/**************************************************************************************************
void DistanceCommand::appendFiles(string temp, string filename) {
	try{
		ofstream output;
		ifstream input;
	
		//open output file in append mode
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
		m->errorOut(e, "DistanceCommand", "appendFiles");
		exit(1);
	}
}
/**************************************************************************************************/


