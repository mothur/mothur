/*
 *  mothurOut.cpp
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"


/******************************************************/
MothurOut* MothurOut::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new MothurOut();
	}
	return _uniqueInstance;
}
/*********************************************************************************************/
set<string> MothurOut::getCurrentTypes()  {
	try {
        
        set<string> types;
        types.insert("fasta");
        types.insert("accnos");
        types.insert("column");
        types.insert("design");
        types.insert("group");
        types.insert("list");
        types.insert("name");
        types.insert("oligos");
        types.insert("order");
        types.insert("ordergroup");
        types.insert("phylip");
        types.insert("qfile");
        types.insert("relabund");
        types.insert("sabund");
        types.insert("rabund");
        types.insert("sff");
        types.insert("shared");
        types.insert("taxonomy");
        types.insert("tree");
        types.insert("flow");
        types.insert("biom");
        types.insert("count");
        types.insert("processors");

		return types;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getCurrentTypes");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::printCurrentFiles()  {
	try {
        
        
		if (accnosfile != "")		{  mothurOut("accnos=" + accnosfile); mothurOutEndLine();			}
		if (columnfile != "")		{  mothurOut("column=" + columnfile); mothurOutEndLine();			}
		if (designfile != "")		{  mothurOut("design=" + designfile); mothurOutEndLine();			}
		if (fastafile != "")		{  mothurOut("fasta=" + fastafile); mothurOutEndLine();				}
		if (groupfile != "")		{  mothurOut("group=" + groupfile); mothurOutEndLine();				}
		if (listfile != "")			{  mothurOut("list=" + listfile); mothurOutEndLine();				}
		if (namefile != "")			{  mothurOut("name=" + namefile); mothurOutEndLine();				}
		if (oligosfile != "")		{  mothurOut("oligos=" + oligosfile); mothurOutEndLine();			}
		if (orderfile != "")		{  mothurOut("order=" + orderfile); mothurOutEndLine();				}
		if (ordergroupfile != "")	{  mothurOut("ordergroup=" + ordergroupfile); mothurOutEndLine();	}
		if (phylipfile != "")		{  mothurOut("phylip=" + phylipfile); mothurOutEndLine();			}
		if (qualfile != "")			{  mothurOut("qfile=" + qualfile); mothurOutEndLine();				}
		if (rabundfile != "")		{  mothurOut("rabund=" + rabundfile); mothurOutEndLine();			}
		if (relabundfile != "")		{  mothurOut("relabund=" + relabundfile); mothurOutEndLine();		}
		if (sabundfile != "")		{  mothurOut("sabund=" + sabundfile); mothurOutEndLine();			}
		if (sfffile != "")			{  mothurOut("sff=" + sfffile); mothurOutEndLine();					}
		if (sharedfile != "")		{  mothurOut("shared=" + sharedfile); mothurOutEndLine();			}
		if (taxonomyfile != "")		{  mothurOut("taxonomy=" + taxonomyfile); mothurOutEndLine();		}
		if (treefile != "")			{  mothurOut("tree=" + treefile); mothurOutEndLine();				}
		if (flowfile != "")			{  mothurOut("flow=" + flowfile); mothurOutEndLine();				}
        if (biomfile != "")			{  mothurOut("biom=" + biomfile); mothurOutEndLine();				}
        if (counttablefile != "")	{  mothurOut("count=" + counttablefile); mothurOutEndLine();	}
		if (processors != "1")		{  mothurOut("processors=" + processors); mothurOutEndLine();		}
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "printCurrentFiles");
		exit(1);
	}
}
/*********************************************************************************************/
bool MothurOut::hasCurrentFiles()  {
	try {
		bool hasCurrent = false;
		
		if (accnosfile != "")		{  return true;			}
		if (columnfile != "")		{  return true;			}
		if (designfile != "")		{  return true;			}
		if (fastafile != "")		{  return true;			}
		if (groupfile != "")		{  return true;			}
		if (listfile != "")			{  return true;			}
		if (namefile != "")			{  return true;			}
		if (oligosfile != "")		{  return true;			}
		if (orderfile != "")		{  return true;			}
		if (ordergroupfile != "")	{  return true;			}
		if (phylipfile != "")		{  return true;			}
		if (qualfile != "")			{  return true;			}
		if (rabundfile != "")		{  return true;			}
		if (relabundfile != "")		{  return true;			}
		if (sabundfile != "")		{  return true;			}
		if (sfffile != "")			{  return true;			}
		if (sharedfile != "")		{  return true;			}
		if (taxonomyfile != "")		{  return true;			}
		if (treefile != "")			{  return true;			}
		if (flowfile != "")			{  return true;			}
        if (biomfile != "")			{  return true;			}
        if (counttablefile != "")	{  return true;			}
		if (processors != "1")		{  return true;			}
		
		return hasCurrent;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "hasCurrentFiles");
		exit(1);
	}
}

/*********************************************************************************************/
void MothurOut::clearCurrentFiles()  {
	try {
		phylipfile = "";
		columnfile = "";
		listfile = "";
		rabundfile = "";
		sabundfile = "";
		namefile = "";
		groupfile = "";
		designfile = "";
		orderfile = "";
		treefile = "";
		sharedfile = "";
		ordergroupfile = "";
		relabundfile = "";
		fastafile = "";
		qualfile = "";
		sfffile = "";
		oligosfile = "";
		accnosfile = "";
		taxonomyfile = "";	
		flowfile = "";
        biomfile = "";
        counttablefile = "";
		processors = "1";
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "clearCurrentFiles");
		exit(1);
	}
}
/***********************************************************************/
string MothurOut::findProgramPath(string programName){
	try { 
		
		string envPath = getenv("PATH");
		string pPath = "";
		
		//delimiting path char
		char delim;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        delim = ':';
#else
        delim = ';';
#endif
		
		//break apart path variable by ':'
		vector<string> dirs;
		splitAtChar(envPath, dirs, delim);
		
        if (debug) { mothurOut("[DEBUG]: dir's in path: \n"); }
        
		//get path related to mothur
		for (int i = 0; i < dirs.size(); i++) {
            
            if (debug) { mothurOut("[DEBUG]: " + dirs[i] + "\n"); }
            
			//to lower so we can find it
			string tempLower = "";
			for (int j = 0; j < dirs[i].length(); j++) {  tempLower += tolower(dirs[i][j]);  }
			
			//is this mothurs path?
			if (tempLower.find(programName) != -1) {  pPath = dirs[i]; break;  }
		}
        
		if (debug) { mothurOut("[DEBUG]: programPath = " + pPath + "\n"); }
        
		if (pPath != "") {
			//add programName so it looks like what argv would look like
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            pPath += "/" + programName;
#else
            pPath += "\\" + programName;
#endif
		}else {
			//okay programName is not in the path, so the folder programName is in must be in the path
			//lets find out which one
			
			//get path related to the program
			for (int i = 0; i < dirs.size(); i++) {
                
                if (debug) { mothurOut("[DEBUG]: looking in " + dirs[i] + " for " + programName + " \n"); }
                
				//is this the programs path?
				ifstream in;
				string tempIn = dirs[i];
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                tempIn += "/" + programName;
#else
                tempIn += "\\" + programName;
#endif
				openInputFile(tempIn, in, "");
				
				//if this file exists
				if (in) { in.close(); pPath = tempIn; if (debug) { mothurOut("[DEBUG]: found it, programPath = " + pPath + "\n"); } break;   }
			}
		}
		
		return pPath;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "findProgramPath");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::setFileName(string filename)  {
	try {
		logFileName = filename;
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		openOutputFile(filename, out);
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setFileName");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::setDefaultPath(string pathname)  {
	try {
	
		//add / to name if needed
		string lastChar = pathname.substr(pathname.length()-1);
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			if (lastChar != "/") { pathname += "/"; }
		#else
			if (lastChar != "\\") { pathname += "\\"; }	
		#endif
		
		defaultPath = pathname;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setDefaultPath");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::setOutputDir(string pathname)  {
	try {
		outputDir = pathname;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setOutputDir");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::closeLog()  {
	try {
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out.close();
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "closeLog");
		exit(1);
	}
}

/*********************************************************************************************/
MothurOut::~MothurOut() {
	try {
		_uniqueInstance = 0;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOut(string output) {
	try {
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out << output;
        logger() << output;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutEndLine() {
	try {
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out << endl;
        logger() << endl;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOut(string output, ofstream& outputFile) {
	try {
		
#ifdef USE_MPI
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if (pid == 0) { //only one process should output to screen
#endif
			
			
			out << output;
			outputFile << output;
            logger() << output;
			
#ifdef USE_MPI
		}
#endif
        
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutEndLine(ofstream& outputFile) {
	try {
#ifdef USE_MPI
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if (pid == 0) { //only one process should output to screen
#endif
			
			out << endl;
			outputFile << endl;
            logger() << endl;
			
#ifdef USE_MPI
		}
#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToLog(string output) {
	try {
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out << output;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutJustToLog");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::errorOut(exception& e, string object, string function) {
	//double vm, rss;
	//mem_usage(vm, rss);
	
	mothurOut("[ERROR]: ");
	mothurOut(toString(e.what()));
	mothurOut(" has occurred in the " + object + " class function " + function + ". Please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
	mothurOutEndLine();
}
/*********************************************************************************************/
//The following was originally from http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c 
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
int MothurOut::mem_usage(double& vm_usage, double& resident_set) {
  #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
  
	   vm_usage     = 0.0;
	   resident_set = 0.0;

	   // 'file' stat seems to give the most reliable results
	   //
	   ifstream stat_stream("/proc/self/stat",ios_base::in);

	   // dummy vars for leading entries in stat that we don't care about
	   //
	   string pid, comm, state, ppid, pgrp, session, tty_nr;
	   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	   string utime, stime, cutime, cstime, priority, nice;
	   string O, itrealvalue, starttime;

	   // the two fields we want
	   //
	   unsigned long vsize;
	   long rss;

	   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
				   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
				   >> utime >> stime >> cutime >> cstime >> priority >> nice
				   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	   vm_usage     = vsize / 1024.0;
	   resident_set = rss * page_size_kb;
	   
	   mothurOut("Memory Usage: vm = " + toString(vm_usage) + " rss = " + toString(resident_set) + "\n");
		return 0;

	#else
/*		//windows memory usage
		// Get the list of process identifiers.
		DWORD aProcesses[1024], cbNeeded, cProcesses;
		
		if ( !EnumProcesses( aProcesses, sizeof(aProcesses), &cbNeeded ) ){ return 1; }

		// Calculate how many process identifiers were returned.
		cProcesses = cbNeeded / sizeof(DWORD);

		// Print the memory usage for each process
		for (int i = 0; i < cProcesses; i++ ) {
			DWORD processID = aProcesses[i];
			
			PROCESS_MEMORY_COUNTERS pmc;

			HANDLE hProcess = OpenProcess((PROCESS_QUERY_INFORMATION | PROCESS_VM_READ), FALSE, processID);

			// Print the process identifier.
			printf( "\nProcess ID: %u\n", processID);
			
			if (NULL != hProcess) {

				if ( GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc)) ) {
					printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
					printf( "\tPeakWorkingSetSize: 0x%08X\n", pmc.PeakWorkingSetSize );
					printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
					printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakPagedPoolUsage );
					printf( "\tQuotaPagedPoolUsage: 0x%08X\n", pmc.QuotaPagedPoolUsage );
					printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakNonPagedPoolUsage );
					printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n", pmc.QuotaNonPagedPoolUsage );
					printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage ); 
					printf( "\tPeakPagefileUsage: 0x%08X\n", pmc.PeakPagefileUsage );
				}
				CloseHandle(hProcess);
			}
		}
*/
			return 0;

	#endif
}


/***********************************************************************/
int MothurOut::openOutputFileAppend(string fileName, ofstream& fileHandle){
	try {
		fileName = getFullPathName(fileName);
		
		fileHandle.open(fileName.c_str(), ios::app);
		if(!fileHandle) {
			mothurOut("[ERROR]: Could not open " + fileName); mothurOutEndLine();
			return 1;
		}
		else {
			return 0;
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "openOutputFileAppend");
		exit(1);
	}
}
/***********************************************************************/
void MothurOut::gobble(istream& f){
	try {
		
		char d;
		while(isspace(d=f.get()))		{ ;}
		if(!f.eof()) { f.putback(d); }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "gobble");
		exit(1);
	}
}
/***********************************************************************/
void MothurOut::gobble(istringstream& f){
	try {
		char d;
		while(isspace(d=f.get()))		{;}
		if(!f.eof()) { f.putback(d); }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "gobble");
		exit(1);
	}
}

/***********************************************************************/

string MothurOut::getline(istringstream& fileHandle) {
	try {
	
		string line = "";
		
		while (!fileHandle.eof())	{
			//get next character
			char c = fileHandle.get(); 
			
			//are you at the end of the line
			if ((c == '\n') || (c == '\r') || (c == '\f')){  break;	}	
			else {		line += c;		}
		}
		
		return line;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getline");
		exit(1);
	}
}
/***********************************************************************/

string MothurOut::getline(ifstream& fileHandle) {
	try {
	
		string line = "";
		
		while (fileHandle)	{
			//get next character
			char c = fileHandle.get(); 
			
			//are you at the end of the line
			if ((c == '\n') || (c == '\r') || (c == '\f') || (c == EOF)){  break;	}	
			else {		line += c;		}
		}
		
		return line;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getline");
		exit(1);
	}
}
/***********************************************************************/

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#ifdef USE_COMPRESSION
inline bool endsWith(string s, const char * suffix){
  size_t suffixLength = strlen(suffix);
  return s.size() >= suffixLength && s.substr(s.size() - suffixLength, suffixLength).compare(suffix) == 0;
}
#endif
#endif

string MothurOut::getRootName(string longName){
	try {
	
		string rootName = longName;

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#ifdef USE_COMPRESSION
    if (endsWith(rootName, ".gz") || endsWith(rootName, ".bz2")) {
      int pos = rootName.find_last_of('.');
      rootName = rootName.substr(0, pos);
      cerr << "shortening " << longName << " to " << rootName << "\n";
    }
#endif
#endif
		if(rootName.find_last_of(".") != rootName.npos){
			int pos = rootName.find_last_of('.')+1;
			rootName = rootName.substr(0, pos);
		}

		return rootName;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getRootName");
		exit(1);
	}
}
/***********************************************************************/

string MothurOut::getSimpleName(string longName){
	try {
		string simpleName = longName;
		
		size_t found;
		found=longName.find_last_of("/\\");

		if(found != longName.npos){
			simpleName = longName.substr(found+1);
		}
		
		return simpleName;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getSimpleName");
		exit(1);
	}
}

/***********************************************************************/

int MothurOut::getRandomIndex(int highest){
	try {
		
		int random = (int) ((float)(highest+1) * (float)(rand()) / ((float)RAND_MAX+1.0));
		
		return random;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getRandomIndex");
		exit(1);
	}	
	
}
/**********************************************************************/

string MothurOut::getPathName(string longName){
	try {
		string rootPathName = longName;
		
		if(longName.find_last_of("/\\") != longName.npos){
			int pos = longName.find_last_of("/\\")+1;
			rootPathName = longName.substr(0, pos);
		}
		
		return rootPathName;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getPathName");
		exit(1);
	}	

}
/***********************************************************************/

bool MothurOut::dirCheck(string& dirName){
	try {
        
        string tag = "";
        #ifdef USE_MPI
            int pid; 
            MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
            tag = toString(pid);
        #endif

        //add / to name if needed
        string lastChar = dirName.substr(dirName.length()-1);
        #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        if (lastChar != "/") { dirName += "/"; }
        #else
        if (lastChar != "\\") { dirName += "\\"; }	
        #endif

        //test to make sure directory exists
        dirName = getFullPathName(dirName);
        string outTemp = dirName + tag + "temp";
        ofstream out;
        out.open(outTemp.c_str(), ios::trunc);
        if(!out) {
            mothurOut(dirName + " directory does not exist or is not writable."); mothurOutEndLine(); 
        }else{
            out.close();
            mothurRemove(outTemp);
            return true;
        }
        
        return false;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "dirCheck");
		exit(1);
	}	
    
}
/***********************************************************************/

string MothurOut::hasPath(string longName){
	try {
		string path = "";
		
		size_t found;
		found=longName.find_last_of("~/\\");

		if(found != longName.npos){
			path = longName.substr(0, found+1);
		}
		
		return path;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "hasPath");
		exit(1);
	}	
}

/***********************************************************************/

string MothurOut::getExtension(string longName){
	try {
		string extension = "";
		
		if(longName.find_last_of('.') != longName.npos){
			int pos = longName.find_last_of('.');
			extension = longName.substr(pos, longName.length());
		}
		
		return extension;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getExtension");
		exit(1);
	}	
}
/***********************************************************************/
bool MothurOut::isBlank(string fileName){
	try {
		
		fileName = getFullPathName(fileName);
		
		ifstream fileHandle;
		fileHandle.open(fileName.c_str());
		if(!fileHandle) {
			mothurOut("[ERROR]: Could not open " + fileName); mothurOutEndLine();
			return false;
		}else {
			//check for blank file
			gobble(fileHandle);
			if (fileHandle.eof()) { fileHandle.close(); return true;  }
			fileHandle.close();
		}
		return false;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "isBlank");
		exit(1);
	}	
}
/***********************************************************************/

string MothurOut::getFullPathName(string fileName){
	try{
	
	string path = hasPath(fileName);
	string newFileName;
	int pos;
	
	if (path == "") { return fileName; } //its a simple name
	else { //we need to complete the pathname
		// ex. ../../../filename 
		// cwd = /user/work/desktop
				
		string cwd;
		//get current working directory 
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)	
			
			if (path.find("~") != -1) { //go to home directory
				string homeDir;
			
				char *homepath = NULL;
				homepath = getenv ("HOME");
				if ( homepath != NULL) { homeDir = homepath; }
				else { homeDir = "";  }

				newFileName = homeDir + fileName.substr(fileName.find("~")+1);
				return newFileName;
			}else { //find path
				if (path.rfind("./") == string::npos) { return fileName; } //already complete name
				else { newFileName = fileName.substr(fileName.rfind("./")+2); } //save the complete part of the name
				
				//char* cwdpath = new char[1024];
				//size_t size;
				//cwdpath=getcwd(cwdpath,size);
				//cwd = cwdpath;
				
				char *cwdpath = NULL;
				cwdpath = getcwd(NULL, 0); // or _getcwd
				if ( cwdpath != NULL) { cwd = cwdpath; }
				else { cwd = "";  }

				
				//rip off first '/'
				string simpleCWD;
				if (cwd.length() > 0) { simpleCWD = cwd.substr(1); }
				
				//break apart the current working directory
				vector<string> dirs;
				while (simpleCWD.find_first_of('/') != string::npos) {
					string dir = simpleCWD.substr(0,simpleCWD.find_first_of('/'));
					simpleCWD = simpleCWD.substr(simpleCWD.find_first_of('/')+1, simpleCWD.length());
					dirs.push_back(dir);
				}
				//get last one              // ex. ../../../filename = /user/work/desktop/filename
				dirs.push_back(simpleCWD);  //ex. dirs[0] = user, dirs[1] = work, dirs[2] = desktop
				
			
				int index = dirs.size()-1;
		
				while((pos = path.rfind("./")) != string::npos) { //while you don't have a complete path
					if (pos == 0) { break;  //you are at the end
					}else if (path[(pos-1)] == '.') { //you want your parent directory ../
						path = path.substr(0, pos-1);
						index--;
						if (index == 0) {  break; }
					}else if (path[(pos-1)] == '/') { //you want the current working dir ./
						path = path.substr(0, pos);
					}else if (pos == 1) { break;  //you are at the end
					}else { mothurOut("cannot resolve path for " +  fileName + "\n"); return fileName; }
				}
			
				for (int i = index; i >= 0; i--) {
					newFileName = dirs[i] +  "/" + newFileName;		
				}
				
				newFileName =  "/" +  newFileName;
				return newFileName;
			}	
		#else
			if (path.find("~") != string::npos) { //go to home directory
				string homeDir = getenv ("HOMEPATH");
				newFileName = homeDir + fileName.substr(fileName.find("~")+1);
				return newFileName;
			}else { //find path
				if (path.rfind(".\\") == string::npos) { return fileName; } //already complete name
				else { newFileName = fileName.substr(fileName.rfind(".\\")+2); } //save the complete part of the name
							
				char *cwdpath = NULL;
				cwdpath = getcwd(NULL, 0); // or _getcwd
				if ( cwdpath != NULL) { cwd = cwdpath; }
				else { cwd = "";  }
				
				//break apart the current working directory
				vector<string> dirs;
				while (cwd.find_first_of('\\') != -1) {
					string dir = cwd.substr(0,cwd.find_first_of('\\'));
					cwd = cwd.substr(cwd.find_first_of('\\')+1, cwd.length());
					dirs.push_back(dir);
		
				}
				//get last one
				dirs.push_back(cwd);  //ex. dirs[0] = user, dirs[1] = work, dirs[2] = desktop
					
				int index = dirs.size()-1;
					
				while((pos = path.rfind(".\\")) != string::npos) { //while you don't have a complete path
					if (pos == 0) { break;  //you are at the end
					}else if (path[(pos-1)] == '.') { //you want your parent directory ../
						path = path.substr(0, pos-1);
						index--;
						if (index == 0) {  break; }
					}else if (path[(pos-1)] == '\\') { //you want the current working dir ./
						path = path.substr(0, pos);
					}else if (pos == 1) { break;  //you are at the end
					}else { mothurOut("cannot resolve path for " +  fileName + "\n"); return fileName; }
				}
			
				for (int i = index; i >= 0; i--) {
					newFileName = dirs[i] +  "\\" + newFileName;		
				}
				
				return newFileName;
			}
			
		#endif
	}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getFullPathName");
		exit(1);
	}	
}
/***********************************************************************/

int MothurOut::openInputFile(string fileName, ifstream& fileHandle, string m){
	try {
			//get full path name
			string completeFileName = getFullPathName(fileName);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#ifdef USE_COMPRESSION
      // check for gzipped or bzipped file
      if (endsWith(completeFileName, ".gz") || endsWith(completeFileName, ".bz2")) {
        string tempName = string(tmpnam(0));
        mkfifo(tempName.c_str(), 0666);
        int fork_result = fork();
        if (fork_result < 0) {
          cerr << "Error forking.\n";
          exit(1);
        } else if (fork_result == 0) {
          string command = (endsWith(completeFileName, ".gz") ? "zcat " : "bzcat ") + completeFileName + string(" > ") + tempName;
          cerr << "Decompressing " << completeFileName << " via temporary named pipe " << tempName << "\n";
          system(command.c_str());
          cerr << "Done decompressing " << completeFileName << "\n";
          mothurRemove(tempName);
          exit(EXIT_SUCCESS);
        } else {
          cerr << "waiting on child process " << fork_result << "\n";
          completeFileName = tempName;
        }
      }
#endif
#endif
			fileHandle.open(completeFileName.c_str());
			if(!fileHandle) {
				//mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
				return 1;
			}else {
				//check for blank file
				gobble(fileHandle);
				return 0;
			}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "openInputFile - no Error");
		exit(1);
	}
}
/***********************************************************************/

int MothurOut::openInputFile(string fileName, ifstream& fileHandle){
	try {

		//get full path name
		string completeFileName = getFullPathName(fileName);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#ifdef USE_COMPRESSION
  // check for gzipped or bzipped file
  if (endsWith(completeFileName, ".gz") || endsWith(completeFileName, ".bz2")) {
    string tempName = string(tmpnam(0));
    mkfifo(tempName.c_str(), 0666);
    int fork_result = fork();
    if (fork_result < 0) {
      cerr << "Error forking.\n";
      exit(1);
    } else if (fork_result == 0) {
      string command = (endsWith(completeFileName, ".gz") ? "zcat " : "bzcat ") + completeFileName + string(" > ") + tempName;
      cerr << "Decompressing " << completeFileName << " via temporary named pipe " << tempName << "\n";
      system(command.c_str());
      cerr << "Done decompressing " << completeFileName << "\n";
      mothurRemove(tempName);
      exit(EXIT_SUCCESS);
    } else {
      cerr << "waiting on child process " << fork_result << "\n";
      completeFileName = tempName;
    }
  }
#endif
#endif

		fileHandle.open(completeFileName.c_str());
		if(!fileHandle) {
			mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
			return 1;
		}
		else {
			//check for blank file
			gobble(fileHandle);
			if (fileHandle.eof()) { mothurOut("[ERROR]: " + completeFileName + " is blank. Please correct."); mothurOutEndLine();  }
			
			return 0;
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "openInputFile");
		exit(1);
	}	
}
/***********************************************************************/

int MothurOut::renameFile(string oldName, string newName){
	try {
        
        if (oldName == newName) { return 0; }
        
		ifstream inTest;
		int exist = openInputFile(newName, inTest, "");
		inTest.close();
		
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
		if (exist == 0) { //you could open it so you want to delete it
			string command = "rm " + newName;
			system(command.c_str());
		}
				
		string command = "mv " + oldName + " " + newName;
		system(command.c_str());
	#else
		mothurRemove(newName);
		int renameOk = rename(oldName.c_str(), newName.c_str());
	#endif
		return 0;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "renameFile");
		exit(1);
	}	
}

/***********************************************************************/

int MothurOut::openOutputFile(string fileName, ofstream& fileHandle){
	try { 
	
		string completeFileName = getFullPathName(fileName);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#ifdef USE_COMPRESSION
    // check for gzipped file
    if (endsWith(completeFileName, ".gz") || endsWith(completeFileName, ".bz2")) {
      string tempName = string(tmpnam(0));
      mkfifo(tempName.c_str(), 0666);
      cerr << "Compressing " << completeFileName << " via temporary named pipe " << tempName << "\n";
      int fork_result = fork();
      if (fork_result < 0) {
        cerr << "Error forking.\n";
        exit(1);
      } else if (fork_result == 0) {
        string command = string(endsWith(completeFileName, ".gz") ?  "gzip" : "bzip2") + " -v > " + completeFileName + string(" < ") + tempName;
        system(command.c_str());
        exit(0);
      } else {
        completeFileName = tempName;
      }
    }
#endif
#endif
		fileHandle.open(completeFileName.c_str(), ios::trunc);
		if(!fileHandle) {
			mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
			return 1;
		}
		else {
			return 0;
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "openOutputFile");
		exit(1);
	}	

}

/**************************************************************************************************/
int MothurOut::appendFiles(string temp, string filename) {
	try{
		ofstream output;
		ifstream input;
	
		//open output file in append mode
		openOutputFileAppend(filename, output);
		int ableToOpen = openInputFile(temp, input, "no error");
		//int ableToOpen = openInputFile(temp, input);
		
		int numLines = 0;
		if (ableToOpen == 0) { //you opened it
            
            char buffer[4096];        
            while (!input.eof()) {
                input.read(buffer, 4096);
                output.write(buffer, input.gcount());
                //count number of lines
                for (int i = 0; i < input.gcount(); i++) {  if (buffer[i] == '\n') {numLines++;} }
            }
			input.close();
		}
		
		output.close();
		
		return numLines;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "appendFiles");
		exit(1);
	}	
}

/**************************************************************************************************/
string MothurOut::sortFile(string distFile, string outputDir){
	try {	
	
		//if (outputDir == "") {  outputDir += hasPath(distFile);  }
		string outfile = getRootName(distFile) + "sorted.dist";

		
		//if you can, use the unix sort since its been optimized for years
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			string command = "sort -n -k +3 " + distFile + " -o " + outfile;
			system(command.c_str());
		#else //you are stuck with my best attempt...
			//windows sort does not have a way to specify a column, only a character in the line
			//since we cannot assume that the distance will always be at the the same character location on each line
			//due to variable sequence name lengths, I chose to force the distance into first position, then sort and then put it back.
		
			//read in file line by file and put distance first
			string tempDistFile = distFile + ".temp";
			ifstream input;
			ofstream output;
			openInputFile(distFile, input);
			openOutputFile(tempDistFile, output);

			string firstName, secondName;
			float dist;
			while (input) {
				input >> firstName >> secondName >> dist;
				output << dist << '\t' << firstName << '\t' << secondName << endl;
				gobble(input);
			}
			input.close();
			output.close();
		
	
			//sort using windows sort
			string tempOutfile = outfile + ".temp";
			string command = "sort " + tempDistFile + " /O " + tempOutfile;
			system(command.c_str());
		
			//read in sorted file and put distance at end again
			ifstream input2;
			openInputFile(tempOutfile, input2);
			openOutputFile(outfile, output);
		
			while (input2) {
				input2 >> dist >> firstName >> secondName;
				output << firstName << '\t' << secondName << '\t' << dist << endl;
				gobble(input2);
			}
			input2.close();
			output.close();
		
			//remove temp files
			mothurRemove(tempDistFile);
			mothurRemove(tempOutfile);
		#endif
		
		return outfile;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "sortFile");
		exit(1);
	}	
}
/**************************************************************************************************/
vector<unsigned long long> MothurOut::setFilePosFasta(string filename, int& num) {
	try {
			vector<unsigned long long> positions;
			ifstream inFASTA;
			//openInputFile(filename, inFASTA);
			inFASTA.open(filename.c_str(), ios::binary);
						
			string input;
			unsigned long long count = 0;
			while(!inFASTA.eof()){
				//input = getline(inFASTA); 
				//cout << input << '\t' << inFASTA.tellg() << endl;
				//if (input.length() != 0) {
				//	if(input[0] == '>'){	unsigned long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	 cout << (pos - input.length() - 1) << endl; }
				//}
				//gobble(inFASTA); //has to be here since windows line endings are 2 characters and mess up the positions
				char c = inFASTA.get(); count++;
				if (c == '>') {
					positions.push_back(count-1);
					//cout << count << endl;
				}
			}
			inFASTA.close();
		
			num = positions.size();
		
			/*FILE * pFile;
			long size;
		
			//get num bytes in file
			pFile = fopen (filename.c_str(),"rb");
			if (pFile==NULL) perror ("Error opening file");
			else{
				fseek (pFile, 0, SEEK_END);
				size=ftell (pFile);
				fclose (pFile);
			}*/
			
			unsigned long long size = positions[(positions.size()-1)];
			ifstream in;
			openInputFile(filename, in);
			
			in.seekg(size);
		
			while(in.get()){
				if(in.eof())		{	break;	}
				else				{	size++;	}
			}
			in.close();
		
			positions.push_back(size);
			positions[0] = 0;
		
			return positions;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setFilePosFasta");
		exit(1);
	}
}
/**************************************************************************************************/
vector<unsigned long long> MothurOut::setFilePosEachLine(string filename, int& num) {
	try {
			filename = getFullPathName(filename);
			
			vector<unsigned long long> positions;
			ifstream in;
			//openInputFile(filename, in);
			in.open(filename.c_str(), ios::binary);
		
			string input;
			unsigned long long count = 0;
			positions.push_back(0);
		
			while(!in.eof()){
				//getline counting reads
				char d = in.get(); count++;
				while ((d != '\n') && (d != '\r') && (d != '\f') && (d != in.eof()))	{
					//get next character
					d = in.get(); 
					count++;
				}
				
				if (!in.eof()) {
					d=in.get(); count++;
					while(isspace(d) && (d != in.eof()))		{ d=in.get(); count++;}
				}
				positions.push_back(count-1);
				//cout << count-1 << endl;
			}
			in.close();
		
			num = positions.size()-1;
		
			FILE * pFile;
			unsigned long long size;
			
			//get num bytes in file
			pFile = fopen (filename.c_str(),"rb");
			if (pFile==NULL) perror ("Error opening file");
			else{
				fseek (pFile, 0, SEEK_END);
				size=ftell (pFile);
				fclose (pFile);
			}
		
			positions[(positions.size()-1)] = size;
		
			return positions;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setFilePosEachLine");
		exit(1);
	}
}
/**************************************************************************************************/

vector<unsigned long long> MothurOut::divideFile(string filename, int& proc) {
	try{
		vector<unsigned long long> filePos;
		filePos.push_back(0);
		
		FILE * pFile;
		unsigned long long size;
		
		filename = getFullPathName(filename);
	
		//get num bytes in file
		pFile = fopen (filename.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
		
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				
		//estimate file breaks
		unsigned long long chunkSize = 0;
		chunkSize = size / proc;

		//file to small to divide by processors
		if (chunkSize == 0)  {  proc = 1;	filePos.push_back(size); return filePos;	}
	
		//for each process seekg to closest file break and search for next '>' char. make that the filebreak
		for (int i = 0; i < proc; i++) {
			unsigned long long spot = (i+1) * chunkSize;
			
			ifstream in;
			openInputFile(filename, in);
			in.seekg(spot);
			
			//look for next '>'
			unsigned long long newSpot = spot;
			while (!in.eof()) {
			   char c = in.get();
				
			   if (c == '>') {   in.putback(c); newSpot = in.tellg(); break;  }
			   else if (int(c) == -1) { break; }
				
			}
		
			//there was not another sequence before the end of the file
			unsigned long long sanityPos = in.tellg();

			if (sanityPos == -1) {	break;  }
			else {  filePos.push_back(newSpot);  }
			
			in.close();
		}
		
		//save end pos
		filePos.push_back(size);
		
		//sanity check filePos
		for (int i = 0; i < (filePos.size()-1); i++) {
			if (filePos[(i+1)] <= filePos[i]) {  filePos.erase(filePos.begin()+(i+1)); i--; }
		}

		proc = (filePos.size() - 1);
#else
		mothurOut("[ERROR]: Windows version should not be calling the divideFile function."); mothurOutEndLine();
		proc=1;
		filePos.push_back(size);
#endif
		return filePos;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "divideFile");
		exit(1);
	}
}
/**************************************************************************************************/
int MothurOut::divideFile(string filename, int& proc, vector<string>& files) {
	try{
		
		vector<unsigned long long> filePos = divideFile(filename, proc);
		
		for (int i = 0; i < (filePos.size()-1); i++) {
			
			//read file chunk
			ifstream in;
			openInputFile(filename, in);
			in.seekg(filePos[i]);
			unsigned long long size = filePos[(i+1)] - filePos[i];
			char* chunk = new char[size];
			in.read(chunk, size);
			in.close();
			
			//open new file
			string fileChunkName = filename + "." + toString(i) + ".tmp";
			ofstream out; 
			openOutputFile(fileChunkName, out);
			
			out << chunk << endl;
			out.close();
			delete[] chunk;
			
			//save name
			files.push_back(fileChunkName);
		}
				
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "divideFile");
		exit(1);
	}
}
/***********************************************************************/

bool MothurOut::isTrue(string f){
	try {
		
		for (int i = 0; i < f.length(); i++) { f[i] = toupper(f[i]); }
		
		if ((f == "TRUE") || (f == "T")) {	return true;	}
		else {	return false;  }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "isTrue");
		exit(1);
	}
}

/***********************************************************************/

float MothurOut::roundDist(float dist, int precision){
	try {
		return int(dist * precision + 0.5)/float(precision);
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "roundDist");
		exit(1);
	}
}
/***********************************************************************/

float MothurOut::ceilDist(float dist, int precision){
	try {
		return int(ceil(dist * precision))/float(precision);
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "ceilDist");
		exit(1);
	}
}
/***********************************************************************/

vector<string> MothurOut::splitWhiteSpace(string& rest, char buffer[], int size){
	try {
        vector<string> pieces;
        
        for (int i = 0; i < size; i++) {
            if (!isspace(buffer[i]))  { rest += buffer[i];  }
            else {
                if (rest != "") { pieces.push_back(rest);  rest = ""; }
                while (i < size) {  //gobble white space
                    if (isspace(buffer[i])) { i++; }
                    else { rest = buffer[i];  break; } //cout << "next piece buffer = " << nextPiece << endl;
                } 
            }
        }
        
        return pieces;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitWhiteSpace");
		exit(1);
	}
}
/***********************************************************************/
vector<string> MothurOut::splitWhiteSpace(string input){
	try {
        vector<string> pieces;
        string rest = "";
        
        for (int i = 0; i < input.length(); i++) {
            if (!isspace(input[i]))  { rest += input[i];  }
            else {
                if (rest != "") { pieces.push_back(rest);  rest = ""; }
                while (i < input.length()) {  //gobble white space
                    if (isspace(input[i])) { i++; }
                    else { rest = input[i];  break; } //cout << "next piece buffer = " << nextPiece << endl;
                } 
            }
        }
        
        if (rest != "") { pieces.push_back(rest); }
        
        return pieces;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitWhiteSpace");
		exit(1);
	}
}
//**********************************************************************************************************************
int MothurOut::readTax(string namefile, map<string, string>& taxMap) {
	try {
		
        //open input file
		ifstream in;
		openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    //are there confidence scores, if so remove them
                    if (secondCol.find_first_of('(') != -1) {  removeConfidences(secondCol);	}
                    taxMap[firstCol] = secondCol;
                    if (debug) {  mothurOut("[DEBUG]: name = '" + firstCol + "' tax = '" + secondCol + "'\n");  }
                    pairDone = false; 
                }
            }
		}
		in.close();
		
		return taxMap.size();

	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readTax");
		exit(1);
	}
}
/**********************************************************************************************************************/
int MothurOut::readNames(string namefile, map<string, string>& nameMap, bool redund) { 
	try {
		
		//open input file
		ifstream in;
		openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = firstCol;  }
                    pairDone = false; 
                }
            }
		}
		in.close();
		
		return nameMap.size();
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
int MothurOut::readNames(string namefile, map<string, string>& nameMap, int flip) { 
	try {
		
		//open input file
		ifstream in;
		openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    nameMap[secondCol] = firstCol;
                    pairDone = false; 
                }
            }
		}
		in.close();
		
		return nameMap.size();
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
int MothurOut::readNames(string namefile, map<string, string>& nameMap, map<string, int>& nameCount) { 
	try {
		nameMap.clear(); nameCount.clear();
		//open input file
		ifstream in;
		openInputFile(namefile, in);
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = firstCol;  }
                    nameCount[firstCol] = theseNames.size();
                    pairDone = false; 
                }
            }
		}
		in.close();
		
		return nameMap.size();
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
int MothurOut::readNames(string namefile, map<string, string>& nameMap) { 
	try {
		
		//open input file
		ifstream in;
		openInputFile(namefile, in);

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
             
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { nameMap[firstCol] = secondCol; pairDone = false; }
            }
		}
		in.close();
		
		return nameMap.size();
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
int MothurOut::readNames(string namefile, map<string, vector<string> >& nameMap) { 
	try {
		
		//open input file
		ifstream in;
		openInputFile(namefile, in);
		
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    vector<string> temp;
                    splitAtComma(secondCol, temp);
                    nameMap[firstCol] = temp;
                    pairDone = false;  
                } 
            }
		}
		in.close();
        
		return nameMap.size();
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
map<string, int> MothurOut::readNames(string namefile) { 
	try {
		
		map<string, int> nameMap;
		
		//open input file
		ifstream in;
		openInputFile(namefile, in);
		
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;  
                } 
            }
		}
        in.close();
		
		return nameMap;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
int MothurOut::readNames(string namefile, vector<seqPriorityNode>& nameVector, map<string, string>& fastamap) { 
	try {
		int error = 0;
		
		//open input file
		ifstream in;
		openInputFile(namefile, in);
		
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        string firstCol, secondCol;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    int num = getNumNames(secondCol);
                    
                    map<string, string>::iterator it = fastamap.find(firstCol);
                    if (it == fastamap.end()) {
                        error = 1;
                        mothurOut("[ERROR]: " + firstCol + " is not in your fastafile, but is in your namesfile, please correct."); mothurOutEndLine();
                    }else {
                        seqPriorityNode temp(num, it->second, firstCol);
                        nameVector.push_back(temp);
                    }
                    
                    pairDone = false;  
                } 
            }
		}
        in.close();
        
		return error;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
//**********************************************************************************************************************
set<string> MothurOut::readAccnos(string accnosfile){
	try {
		set<string> names;
		ifstream in;
		openInputFile(accnosfile, in);
		string name;
		
        string rest = "";
        char buffer[4096];
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {  names.insert(pieces[i]);  }
        }
		in.close();	
		
		return names;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readAccnos");
		exit(1);
	}
}
//**********************************************************************************************************************
int MothurOut::readAccnos(string accnosfile, vector<string>& names){
	try {
        names.clear();
		ifstream in;
		openInputFile(accnosfile, in);
		string name;
		
        string rest = "";
        char buffer[4096];
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {  names.push_back(pieces[i]);  }
        }
		in.close();	
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readAccnos");
		exit(1);
	}
}
/***********************************************************************/

int MothurOut::getNumNames(string names){
	try {
		int count = 0;
		
		if(names != ""){
			count = 1;
			for(int i=0;i<names.size();i++){
				if(names[i] == ','){
					count++;
				}
			}
		}
		
		return count;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getNumNames");
		exit(1);
	}
}
/***********************************************************************/

int MothurOut::getNumChar(string line, char c){
	try {
		int count = 0;
		
		if(line != ""){
			for(int i=0;i<line.size();i++){
				if(line[i] == c){
					count++;
				}
			}
		}
		
		return count;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getNumChar");
		exit(1);
	}
}
/***********************************************************************/
int MothurOut::mothurRemove(string filename){
	try {
		filename = getFullPathName(filename);
		int error = remove(filename.c_str());
		//if (error != 0) { 
		//	if (errno != ENOENT) { //ENOENT == file does not exist
		//		string message = "Error deleting file " + filename;
		//		perror(message.c_str()); 
		//	}
		//}
		return error;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "mothurRemove");
		exit(1);
	}
}
/***********************************************************************/
bool MothurOut::mothurConvert(string item, int& num){
	try {
		bool error = false;
		
		if (isNumeric1(item)) {
			convert(item, num);
		}else {
			num = 0;
			error = true;
			mothurOut("[ERROR]: cannot convert " + item + " to an integer."); mothurOutEndLine();
			commandInputsConvertError = true;
		}
		
		return error;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "mothurConvert");
		exit(1);
	}
}
/***********************************************************************/
bool MothurOut::mothurConvert(string item, intDist& num){
	try {
		bool error = false;
		
		if (isNumeric1(item)) {
			convert(item, num);
		}else {
			num = 0;
			error = true;
			mothurOut("[ERROR]: cannot convert " + item + " to an integer."); mothurOutEndLine();
			commandInputsConvertError = true;
		}
		
		return error;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "mothurConvert");
		exit(1);
	}
}

/***********************************************************************/
bool MothurOut::isNumeric1(string stringToCheck){
	try {
		bool numeric = false;
		
		if(stringToCheck.find_first_not_of("0123456789.-") == string::npos) { numeric = true; }
			
		return numeric;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "isNumeric1");
		exit(1);
	}
	
}
/***********************************************************************/
bool MothurOut::mothurConvert(string item, float& num){
	try {
		bool error = false;
		
		if (isNumeric1(item)) {
			convert(item, num);
		}else {
			num = 0;
			error = true;
			mothurOut("[ERROR]: cannot convert " + item + " to a float."); mothurOutEndLine();
			commandInputsConvertError = true;
		}
		
		return error;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "mothurConvert");
		exit(1);
	}
}
/***********************************************************************/
bool MothurOut::mothurConvert(string item, double& num){
	try {
		bool error = false;
		
		if (isNumeric1(item)) {
			convert(item, num);
		}else {
			num = 0;
			error = true;
			mothurOut("[ERROR]: cannot convert " + item + " to a double."); mothurOutEndLine();
			commandInputsConvertError = true;
		}
		
		return error;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "mothurConvert");
		exit(1);
	}
}
/**************************************************************************************************/

vector<vector<double> > MothurOut::binomial(int maxOrder){
	try {
	vector<vector<double> > binomial(maxOrder+1);
	
    for(int i=0;i<=maxOrder;i++){
		binomial[i].resize(maxOrder+1);
		binomial[i][0]=1;
		binomial[0][i]=0;
    }
    binomial[0][0]=1;
	
    binomial[1][0]=1;
    binomial[1][1]=1;
	
    for(int i=2;i<=maxOrder;i++){
		binomial[1][i]=0;
    }
	
    for(int i=2;i<=maxOrder;i++){
		for(int j=1;j<=maxOrder;j++){
			if(i==j){	binomial[i][j]=1;									}
			if(j>i)	{	binomial[i][j]=0;									}
			else	{	binomial[i][j]=binomial[i-1][j-1]+binomial[i-1][j];	}
		}
    }
	
	return binomial;
	
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "binomial");
		exit(1);
	}
}
/**************************************************************************************************/
unsigned int MothurOut::fromBase36(string base36){
	try {
		unsigned int num = 0;
		
		map<char, int> converts;
		converts['A'] = 0;
		converts['a'] = 0;
		converts['B'] = 1;
		converts['b'] = 1;
		converts['C'] = 2;
		converts['c'] = 2;
		converts['D'] = 3;
		converts['d'] = 3;
		converts['E'] = 4;
		converts['e'] = 4;
		converts['F'] = 5;
		converts['f'] = 5;
		converts['G'] = 6;
		converts['g'] = 6;
		converts['H'] = 7;
		converts['h'] = 7;
		converts['I'] = 8;
		converts['i'] = 8;
		converts['J'] = 9;
		converts['j'] = 9;
		converts['K'] = 10;
		converts['k'] = 10;
		converts['L'] = 11;
		converts['l'] = 11;
		converts['M'] = 12;
		converts['m'] = 12;
		converts['N'] = 13;
		converts['n'] = 13;
		converts['O'] = 14;
		converts['o'] = 14;
		converts['P'] = 15;
		converts['p'] = 15;
		converts['Q'] = 16;
		converts['q'] = 16;
		converts['R'] = 17;
		converts['r'] = 17;
		converts['S'] = 18;
		converts['s'] = 18;
		converts['T'] = 19;
		converts['t'] = 19;
		converts['U'] = 20;
		converts['u'] = 20;
		converts['V'] = 21;
		converts['v'] = 21;
		converts['W'] = 22;
		converts['w'] = 22;
		converts['X'] = 23;
		converts['x'] = 23;
		converts['Y'] = 24;
		converts['y'] = 24;
		converts['Z'] = 25;
		converts['z'] = 25;
		converts['0'] = 26;
		converts['1'] = 27;
		converts['2'] = 28;
		converts['3'] = 29;
		converts['4'] = 30;
		converts['5'] = 31;
		converts['6'] = 32;
		converts['7'] = 33;
		converts['8'] = 34;
		converts['9'] = 35;		
		
		int i = 0;
		while (i < base36.length()) {
			char c = base36[i];
			num = 36 * num + converts[c];
			i++;
		}
		
		return num;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "fromBase36");
		exit(1);
	}
}
/***********************************************************************/

int MothurOut::factorial(int num){
	try {
		int total = 1;
		
		for (int i = 1; i <= num; i++) {
			total *= i;
		}
		
		return total;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "factorial");
		exit(1);
	}
}
/***********************************************************************/

int MothurOut::getNumSeqs(ifstream& file){
	try {
		int numSeqs = count(istreambuf_iterator<char>(file),istreambuf_iterator<char>(), '>');
		file.seekg(0);
		return numSeqs;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getNumSeqs");
		exit(1);
	}	
}
/***********************************************************************/
void MothurOut::getNumSeqs(ifstream& file, int& numSeqs){
	try {
		string input;
		numSeqs = 0;
		while(!file.eof()){
			input = getline(file);
			if (input.length() != 0) {
				if(input[0] == '>'){ numSeqs++;	}
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getNumSeqs");
		exit(1);
	}	
}
/***********************************************************************/

//This function parses the estimator options and puts them in a vector
void MothurOut::splitAtChar(string& estim, vector<string>& container, char symbol) {
	try {
		string individual = "";
		int estimLength = estim.size();
		for(int i=0;i<estimLength;i++){
			if(estim[i] == symbol){
				container.push_back(individual);
				individual = "";				
			}
			else{
				individual += estim[i];
			}
		}
		container.push_back(individual);

	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtChar");
		exit(1);
	}	
}

/***********************************************************************/

//This function parses the estimator options and puts them in a vector
void MothurOut::splitAtDash(string& estim, vector<string>& container) {
	try {
		string individual = "";
		int estimLength = estim.size();
		for(int i=0;i<estimLength;i++){
			if(estim[i] == '-'){
				container.push_back(individual);
				individual = "";				
			}
			else{
				individual += estim[i];
			}
		}
		container.push_back(individual);

	
	/*	string individual;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				container.push_back(individual);
			}
		}
		//get last one
		container.push_back(estim); */
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtDash");
		exit(1);
	}	
}

/***********************************************************************/
//This function parses the label options and puts them in a set
void MothurOut::splitAtDash(string& estim, set<string>& container) {
	try {
		string individual = "";
		int estimLength = estim.size();
		for(int i=0;i<estimLength;i++){
			if(estim[i] == '-'){
				container.insert(individual);
				individual = "";				
			}
			else{
				individual += estim[i];
			}
		}
		container.insert(individual);

	//	string individual;
		
	//	while (estim.find_first_of('-') != -1) {
	//		individual = estim.substr(0,estim.find_first_of('-'));
	//		if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
	//			estim = estim.substr(estim.find_first_of('-')+1, estim.length());
	//			container.insert(individual);
	//		}
	//	}
		//get last one
	//	container.insert(estim);
	
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtDash");
		exit(1);
	}	
}
/***********************************************************************/
//This function parses the line options and puts them in a set
void MothurOut::splitAtDash(string& estim, set<int>& container) {
	try {
		string individual;
		int lineNum;
		
		while (estim.find_first_of('-') != -1) {
			individual = estim.substr(0,estim.find_first_of('-'));
			if ((estim.find_first_of('-')+1) <= estim.length()) { //checks to make sure you don't have dash at end of string
				estim = estim.substr(estim.find_first_of('-')+1, estim.length());
				convert(individual, lineNum); //convert the string to int
				container.insert(lineNum);
			}
		}
		//get last one
		convert(estim, lineNum); //convert the string to int
		container.insert(lineNum);
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtDash");
		exit(1);
	}	
}
/***********************************************************************/
string MothurOut::makeList(vector<string>& names) {
	try {
		string list = "";
        
        if (names.size() == 0) { return list; }
		
        for (int i = 0; i < names.size()-1; i++) { list += names[i] + ",";  }
        
        //get last name
        list += names[names.size()-1];
        
        return list;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "makeList");
		exit(1);
	}	
}

/***********************************************************************/
//This function parses the a string and puts peices in a vector
void MothurOut::splitAtComma(string& estim, vector<string>& container) {
	try {
		string individual = "";
		int estimLength = estim.size();
		for(int i=0;i<estimLength;i++){
			if(estim[i] == ','){
				container.push_back(individual);
				individual = "";				
			}
			else{
				individual += estim[i];
			}
		}
		container.push_back(individual);
		
		
		
		
//		string individual;
//		
//		while (estim.find_first_of(',') != -1) {
//			individual = estim.substr(0,estim.find_first_of(','));
//			if ((estim.find_first_of(',')+1) <= estim.length()) { //checks to make sure you don't have comma at end of string
//				estim = estim.substr(estim.find_first_of(',')+1, estim.length());
//				container.push_back(individual);
//			}
//		}
//		//get last one
//		container.push_back(estim);
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtComma");
		exit(1);
	}	
}
/***********************************************************************/
//This function splits up the various option parameters
void MothurOut::splitAtChar(string& prefix, string& suffix, char c){
	try {
		prefix = suffix.substr(0,suffix.find_first_of(c));
		if ((suffix.find_first_of(c)+2) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
			suffix = suffix.substr(suffix.find_first_of(c)+1, suffix.length());
			string space = " ";
			while(suffix.at(0) == ' ')
				suffix = suffix.substr(1, suffix.length());
		}
        
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtComma");
		exit(1);
	}	
}

/***********************************************************************/

//This function splits up the various option parameters
void MothurOut::splitAtComma(string& prefix, string& suffix){
	try {
		prefix = suffix.substr(0,suffix.find_first_of(','));
		if ((suffix.find_first_of(',')+2) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
			suffix = suffix.substr(suffix.find_first_of(',')+1, suffix.length());
			string space = " ";
			while(suffix.at(0) == ' ')
				suffix = suffix.substr(1, suffix.length());
		}

	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtComma");
		exit(1);
	}	
}
/***********************************************************************/

//This function separates the key value from the option value i.e. dist=96_...
void MothurOut::splitAtEquals(string& key, string& value){		
	try {
		if(value.find_first_of('=') != -1){
			key = value.substr(0,value.find_first_of('='));
			if ((value.find_first_of('=')+1) <= value.length()) {
				value = value.substr(value.find_first_of('=')+1, value.length());
			}
		}else{
			key = value;
			value = 1;
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtEquals");
		exit(1);
	}	
}

/**************************************************************************************************/

bool MothurOut::inUsersGroups(string groupname, vector<string> Groups) {
	try {
		for (int i = 0; i < Groups.size(); i++) {
			if (groupname == Groups[i]) { return true; }
		}
		return false;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "inUsersGroups");
		exit(1);
	}	
}
/**************************************************************************************************/
//returns true if any of the strings in first vector are in second vector
bool MothurOut::inUsersGroups(vector<string> groupnames, vector<string> Groups) {
	try {
		
		for (int i = 0; i < groupnames.size(); i++) {
			if (inUsersGroups(groupnames[i], Groups)) { return true; }
		}
		return false;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "inUsersGroups");
		exit(1);
	}	
}
/***********************************************************************/
//this function determines if the user has given us labels that are smaller than the given label.
//if so then it returns true so that the calling function can run the previous valid distance.
//it's a "smart" distance function.  It also checks for invalid labels.
bool MothurOut::anyLabelsToProcess(string label, set<string>& userLabels, string errorOff) {
	try {
		
		set<string>::iterator it;
		vector<float> orderFloat;
		map<string, float> userMap;  //the conversion process removes trailing 0's which we need to put back
		map<string, float>::iterator it2;
		float labelFloat;
		bool smaller = false;
		
		//unique is the smallest line
		if (label == "unique") {  return false;  }
		else { 
			if (convertTestFloat(label, labelFloat)) {
				convert(label, labelFloat); 
			}else { //cant convert 
				return false;
			}
		}
		
		//go through users set and make them floats
		for(it = userLabels.begin(); it != userLabels.end();) {
			
			float temp;
			if ((*it != "unique") && (convertTestFloat(*it, temp) == true)){
				convert(*it, temp);
				orderFloat.push_back(temp);
				userMap[*it] = temp;
				it++;
			}else if (*it == "unique") { 
				orderFloat.push_back(-1.0);
				userMap["unique"] = -1.0;
				it++;
			}else {
				if (errorOff == "") {  mothurOut(*it + " is not a valid label."); mothurOutEndLine();  }
				userLabels.erase(it++); 
			}
		}
		
		//sort order
		sort(orderFloat.begin(), orderFloat.end());
		
		/*************************************************/
		//is this label bigger than any of the users labels
		/*************************************************/
				
		//loop through order until you find a label greater than label
		for (int i = 0; i < orderFloat.size(); i++) {
			if (orderFloat[i] < labelFloat) {
				smaller = true;
				if (orderFloat[i] == -1) { 
					if (errorOff == "") { mothurOut("Your file does not include the label unique."); mothurOutEndLine(); }
					userLabels.erase("unique");
				}
				else {  
					if (errorOff == "") { mothurOut("Your file does not include the label "); mothurOutEndLine(); }
					string s = "";
					for (it2 = userMap.begin(); it2!= userMap.end(); it2++) {  
						if (it2->second == orderFloat[i]) {  
							s = it2->first;  
							//remove small labels
							userLabels.erase(s);
							break;
						}
					}
					if (errorOff == "") {mothurOut( s +  ". I will use the next smallest distance. "); mothurOutEndLine(); }
				}
			//since they are sorted once you find a bigger one stop looking
			}else { break; }
		}
		
		return smaller;
						
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "anyLabelsToProcess");
		exit(1);
	}	
}

/**************************************************************************************************/
bool MothurOut::checkReleaseVersion(ifstream& file, string version) {
	try {
		
		bool good = true;
		
		string line = getline(file);  

		//before we added this check
		if (line[0] != '#') {  good = false;  }
		else {
			//rip off #
			line = line.substr(1);
			
			vector<string> versionVector;
			splitAtChar(version, versionVector, '.');
			
			//check file version
			vector<string> linesVector;
			splitAtChar(line, linesVector, '.');
			
			if (versionVector.size() != linesVector.size()) { good = false; }
			else {
				for (int j = 0; j < versionVector.size(); j++) {
					int num1, num2;
					convert(versionVector[j], num1);
					convert(linesVector[j], num2);
					
					//if mothurs version is newer than this files version, then we want to remake it
					if (num1 > num2) {  good = false; break;  }
				}
			}
			
		}
		
		if (!good) {  file.close();  }
		else { file.seekg(0);  }
		
		return good;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "checkReleaseVersion");		
		exit(1);
	}
}
/**************************************************************************************************/
bool MothurOut::isContainingOnlyDigits(string input) {
	try{
		
		//are you a digit in ascii code
		for (int i = 0;i < input.length(); i++){
			if( input[i]>47 && input[i]<58){}
			else { return false; }
		}
		
		return true;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "isContainingOnlyDigits");		
		exit(1);
	}
}
/**************************************************************************************************/
int MothurOut::removeConfidences(string& tax) {
	try {
		
		string taxon;
		string newTax = "";
		
		while (tax.find_first_of(';') != -1) {
			
			if (control_pressed) { return 0; }
			
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'));
	
			int pos = taxon.find_last_of('(');
			if (pos != -1) {
				//is it a number?
				int pos2 = taxon.find_last_of(')');
				if (pos2 != -1) {
					string confidenceScore = taxon.substr(pos+1, (pos2-(pos+1)));
					if (isNumeric1(confidenceScore)) {
						taxon = taxon.substr(0, pos); //rip off confidence 
					}
				}
			}
			taxon += ";";
			
			tax = tax.substr(tax.find_first_of(';')+1, tax.length());
			newTax += taxon;
		}
		
		tax = newTax;
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "removeConfidences");
		exit(1);
	}
}
/**************************************************************************************************/





