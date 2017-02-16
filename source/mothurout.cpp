/*
 *  mothurOut.cpp
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"

//needed for testing project
//MothurOut* MothurOut::_uniqueInstance;

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
        types.insert("summary");
        types.insert("file");
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
        types.insert("constaxonomy");
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
void MothurOut::printCurrentFiles(string filename)  {
	try {
        
        if (filename != "") {
            ofstream out;
            openOutputFile(filename, out);
            
            if (accnosfile != "")		{  mothurOut("accnos=" + accnosfile, out); mothurOutEndLine(out);           }
            if (columnfile != "")		{  mothurOut("column=" + columnfile, out); mothurOutEndLine(out);			}
            if (designfile != "")		{  mothurOut("design=" + designfile, out); mothurOutEndLine(out);			}
            if (fastafile != "")		{  mothurOut("fasta=" + fastafile, out); mothurOutEndLine(out);				}
            if (groupfile != "")		{  mothurOut("group=" + groupfile, out); mothurOutEndLine(out);				}
            if (listfile != "")			{  mothurOut("list=" + listfile, out); mothurOutEndLine(out);				}
            if (namefile != "")			{  mothurOut("name=" + namefile, out); mothurOutEndLine(out);				}
            if (oligosfile != "")		{  mothurOut("oligos=" + oligosfile, out); mothurOutEndLine(out);			}
            if (orderfile != "")		{  mothurOut("order=" + orderfile, out); mothurOutEndLine(out);				}
            if (ordergroupfile != "")	{  mothurOut("ordergroup=" + ordergroupfile, out); mothurOutEndLine(out);	}
            if (phylipfile != "")		{  mothurOut("phylip=" + phylipfile, out); mothurOutEndLine(out);			}
            if (qualfile != "")			{  mothurOut("qfile=" + qualfile, out); mothurOutEndLine(out);				}
            if (rabundfile != "")		{  mothurOut("rabund=" + rabundfile, out); mothurOutEndLine(out);			}
            if (relabundfile != "")		{  mothurOut("relabund=" + relabundfile, out); mothurOutEndLine(out);		}
            if (sabundfile != "")		{  mothurOut("sabund=" + sabundfile, out); mothurOutEndLine(out);			}
            if (sfffile != "")			{  mothurOut("sff=" + sfffile, out); mothurOutEndLine(out);					}
            if (sharedfile != "")		{  mothurOut("shared=" + sharedfile, out); mothurOutEndLine(out);			}
            if (taxonomyfile != "")		{  mothurOut("taxonomy=" + taxonomyfile, out); mothurOutEndLine(out);		}
            if (constaxonomyfile != "")	{  mothurOut("constaxonomy=" + constaxonomyfile, out); mothurOutEndLine(out);}
            if (treefile != "")			{  mothurOut("tree=" + treefile, out); mothurOutEndLine(out);				}
            if (flowfile != "")			{  mothurOut("flow=" + flowfile, out); mothurOutEndLine(out);				}
            if (biomfile != "")			{  mothurOut("biom=" + biomfile, out); mothurOutEndLine(out);				}
            if (counttablefile != "")	{  mothurOut("count=" + counttablefile, out); mothurOutEndLine(out);        }
            if (processors != "1")		{  mothurOut("processors=" + processors, out); mothurOutEndLine(out);		}
            if (summaryfile != "")		{  mothurOut("summary=" + summaryfile, out); mothurOutEndLine(out);         }
            if (filefile != "")         {  mothurOut("file=" + filefile, out); mothurOutEndLine(out);               }
            
            out.close();
            
        }else {
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
            if (constaxonomyfile != "")	{  mothurOut("constaxonomy=" + constaxonomyfile); mothurOutEndLine();}
            if (treefile != "")			{  mothurOut("tree=" + treefile); mothurOutEndLine();				}
            if (flowfile != "")			{  mothurOut("flow=" + flowfile); mothurOutEndLine();				}
            if (biomfile != "")			{  mothurOut("biom=" + biomfile); mothurOutEndLine();				}
            if (counttablefile != "")	{  mothurOut("count=" + counttablefile); mothurOutEndLine();        }
            if (processors != "1")		{  mothurOut("processors=" + processors); mothurOutEndLine();		}
            if (summaryfile != "")		{  mothurOut("summary=" + summaryfile); mothurOutEndLine();         }
            if (filefile != "")         {  mothurOut("file=" + filefile); mothurOutEndLine();               }
        }
		
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
        if (constaxonomyfile != "")	{  return true;			}
		if (treefile != "")			{  return true;			}
		if (flowfile != "")			{  return true;			}
        if (biomfile != "")			{  return true;			}
        if (counttablefile != "")	{  return true;			}
        if (summaryfile != "")      {  return true;			}
        if (filefile != "")         {  return true;			}
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
        filefile = "";
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
        summaryfile = "";
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
        string pPath = "";
        
        //look in ./
        //is this the programs path?
        ifstream in5;
        string tempIn = ".";
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        tempIn += "/" + programName;
#else
        tempIn += "\\" + programName;
#endif
        openInputFile(tempIn, in5, "");
        
        //if this file exists
        if (in5) { in5.close(); pPath = getFullPathName(tempIn); if (debug) { mothurOut("[DEBUG]: found it, programPath = " + pPath + "\n"); } return pPath;   }
		
		string envPath = getenv("PATH");
		
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
		openOutputFile(filename, out);
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
		
		defaultPath = getFullPathName(pathname);
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setDefaultPath");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::setBlastPath(string pathname)  {
    try {
        
        //add / to name if needed
        string lastChar = pathname.substr(pathname.length()-1);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        if (lastChar != "/") { pathname += "/"; }
#else
        if (lastChar != "\\") { pathname += "\\"; }
#endif
        
        blastPath = getFullPathName(pathname);
        
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
        if (numErrors != 0) {
            out << "\n\n************************************************************\n";
            out << "************************************************************\n";
            out << "************************************************************\n";
            out << "Detected " + toString(numErrors) + " [ERROR] messages, please review.\n";
            out << "************************************************************\n";
            out << "************************************************************\n";
            out << "************************************************************\n";
            
            logger() << "\n\n************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "Detected " + toString(numErrors) + " [ERROR] messages, please review.\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
        }
        
        if (numWarnings != 0) {
            out << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            out << "Detected " + toString(numWarnings) + " [WARNING] messages, please review.\n";
            out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            
            logger() << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "Detected " + toString(numWarnings) + " [WARNING] messages, please review.\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            logger() << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
        }
        
		out.close();
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
        if (output.find("[ERROR]") != string::npos) { numErrors++; }
        
        if (output.find("[WARNING]") != string::npos) { numWarnings++; }

        
        if (!quietMode) {
            out << output;
            logger() << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                out << output;
                logger() << output;
            }
        }
        
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToScreen(string output) {
	try {
		if (output.find("[ERROR]") != string::npos) { numErrors++; }
        
        if (output.find("[WARNING]") != string::npos) { numWarnings++; }
        
        if (!quietMode) {
            logger() << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                logger() << output;
            }
        }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutEndLine() {
	try {
		if (!quietMode) {
            out << endl;
            logger() << endl;
        }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOut(string output, ofstream& outputFile) {
	try {
        if (output.find("[ERROR]") != string::npos) { numErrors++; }
        
        if (output.find("[WARNING]") != string::npos) { numWarnings++; }
        
        if (!quietMode) {
            out << output;
            outputFile << output;
            logger() << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                out << output;
                outputFile << output;
                logger() << output;
            }
            
        }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutEndLine(ofstream& outputFile) {
	try {
        if (!quietMode) {
            out << endl;
            outputFile << endl;
            logger() << endl;
        }
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToLog(string output) {
	try {
        if (output.find("[ERROR]") != string::npos) { numErrors++; }

        if (output.find("[WARNING]") != string::npos) { numWarnings++; }
        
        if (!quietMode) {
            out << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                out << output;
            }
        }

	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutJustToLog");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::errorOut(exception& e, string object, string function) {
    numErrors++; 
	
    string errorType = toString(e.what());
    
    int pos = errorType.find("bad_alloc");
    mothurOut("[ERROR]: ");
    mothurOut(errorType);
    
    unsigned long long ramUsed, total;
    ramUsed = getRAMUsed(); total = getTotalRAM();
    mothurOut("RAM used: " + toString(ramUsed/(double)GIG) + "Gigabytes . Total Ram: " + toString(total/(double)GIG) + "Gigabytes.\n\n");
    
    if (pos == string::npos) { //not bad_alloc
        mothurOut(" has occurred in the " + object + " class function " + function + ". Please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
        mothurOutEndLine();
    }else { //bad alloc
        if (object == "cluster"){
            mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory.  There are two common causes for this, file size and format.\n\nFile Size:\nThe cluster command loads your distance matrix into RAM, and your distance file is most likely too large to fit in RAM. There are two options to help with this. The first is to use a cutoff. By using a cutoff mothur will only load distances that are below the cutoff. If that is still not enough, there is a command called cluster.split, http://www.mothur.org/wiki/cluster.split which divides the distance matrix, and clusters the smaller pieces separately. You may also be able to reduce the size of the original distance matrix by using the commands outlined in the Schloss SOP, http://www.mothur.org/wiki/Schloss_SOP. \n\nWrong Format:\nThis error can be caused by trying to read a column formatted distance matrix using the phylip parameter. By default, the dist.seqs command generates a column formatted distance matrix. To make a phylip formatted matrix set the dist.seqs command parameter output to lt.  \n\nIf you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
        }else if (object == "shhh.flows"){
                mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory. The shhh.flows command is very memory intensive. This error is most commonly caused by trying to process a dataset too large, using multiple processors, or failing to run trim.flows before shhh.flows. If you are running our 32bit version, your memory usage is limited to 4G.  If you have more than 4G of RAM and are running a 64bit OS, using our 64bit version may resolve your issue.  If you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required. Running trim.flows with an oligos file, and then shhh.flows with the file option may also resolve the issue. If for some reason you are unable to run shhh.flows with your data, a good alternative is to use the trim.seqs command using a 50-bp sliding window and to trim the sequence when the average quality score over that window drops below 35. Our results suggest that the sequencing error rates by this method are very good, but not quite as good as by shhh.flows and that the resulting sequences tend to be a bit shorter. If you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry. ");
        }else {
            mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory.  This is most commonly caused by trying to process a dataset too large, using multiple processors, or a file format issue. If you are running our 32bit version, your memory usage is limited to 4G.  If you have more than 4G of RAM and are running a 64bit OS, using our 64bit version may resolve your issue.  If you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required. Also, you may be able to reduce the size of your dataset by using the commands outlined in the Schloss SOP, http://www.mothur.org/wiki/Schloss_SOP. If you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
        }
    }
}
/*********************************************************************************************/
unsigned long long MothurOut::getRAMUsed() {
    try {
       
#if defined (__APPLE__) || (__MACH__)
    /* Mac: ru_maxrss gives the size in bytes */
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, & r_usage);
    
    return r_usage.ru_maxrss;
#elif (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
    /* Linux: ru_maxrss gives the size in kilobytes  */
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, & r_usage);
    return r_usage.ru_maxrss * 1024;
#else
        MEMORYSTATUSEX status;
        status.dwLength = sizeof(status);
        GlobalMemoryStatusEx(&status);
        return (size_t)(status.ullTotalPhys - status.ullAvailPhys);
#endif
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "getMemoryUsed");
        exit(1);
    }
}
/*********************************************************************************************/

unsigned long long MothurOut::getTotalRAM() {
    try {
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
    #if defined _SC_PHYS_PAGES && defined _SC_PAGESIZE
        /* This works on linux-gnu, solaris2 and cygwin.  */
        double pages = sysconf (_SC_PHYS_PAGES);
        double pagesize = sysconf (_SC_PAGESIZE);
        if (0 <= pages && 0 <= pagesize)
        return pages * pagesize;
    #else
        mothurOut("[WARNING]: Cannot determine amount of RAM");
    #endif
        
#elif defined (_WIN32)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return (size_t)status.ullTotalPhys;
#else
    struct sysinfo si;
    if (sysinfo(&si))
        mothurOut("[WARNING]: Cannot determine amount of RAM");
    return si.totalram * si.mem_unit;
    
#endif
        return 0;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "getTotalRAM");
        exit(1);
    }
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
int MothurOut::openOutputFileBinaryAppend(string fileName, ofstream& fileHandle){
	try {
		fileName = getFullPathName(fileName);
		
		fileHandle.open(fileName.c_str(), ios::app | ios::binary);
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
void MothurOut::zapGremlins(istream& f){
    try {
        
        char d;
        while('\0'==(d=f.get()))		{ ;}
        if(!f.eof()) { f.putback(d); }
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "zapGremlins");
        exit(1);
    }
}
/***********************************************************************/
void MothurOut::zapGremlins(istringstream& f){
    try {
        char d;
        while('\0'==(d=f.get()))		{ ;}
        if(!f.eof()) { f.putback(d); }
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "zapGremlins");
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
string MothurOut::getRootName(string longName){
	try {
	
		string rootName = longName;

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
        if (highest == 0) { return 0; }
        
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
        
        if (dirName == "") { return false; }
        
        string tag = "";

        //add / to name if needed
        string lastChar = dirName.substr(dirName.length()-1);
        #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        if (lastChar != "/") { dirName += "/"; }
        #else
        if (lastChar != "\\") { dirName += "\\"; }	
        #endif

        //test to make sure directory exists
        dirName = getFullPathName(dirName);
        string outTemp = dirName + tag + "temp"+ toString(time(NULL));
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

bool MothurOut::dirCheck(string& dirName, string noError){
    try {
        
        if (dirName == "") { return false; }
        
        string tag = "";
        
        //add / to name if needed
        string lastChar = dirName.substr(dirName.length()-1);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        if (lastChar != "/") { dirName += "/"; }
#else
        if (lastChar != "\\") { dirName += "\\"; }
#endif
        
        //test to make sure directory exists
        dirName = getFullPathName(dirName);
        string outTemp = dirName + tag + "temp"+ toString(time(NULL));
        ofstream out;
        out.open(outTemp.c_str(), ios::trunc);
        if(!out) {
            //mothurOut(dirName + " directory does not exist or is not writable."); mothurOutEndLine();
        }else{
            out.close();
            mothurRemove(outTemp);
            return true;
        }
        
        return false;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "dirCheck - noError");
        exit(1);
    }
    
}
/***********************************************************************/
//returns true it exits or if we can make it
bool MothurOut::mkDir(string& dirName){
    try {
        
        bool dirExist = dirCheck(dirName, "noError");
        
        if (dirExist) { return true; }
        
        string tag = "";
        
        string makeDirectoryCommand = "";

        makeDirectoryCommand += "mkdir -p \"" + dirName + "\"";
        system(makeDirectoryCommand.c_str());
        if (dirCheck(dirName)) { return true; }
        
        return false;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "mkDir");
        exit(1);
    }
    
}
//**********************************************************************************************************************

map<string, vector<string> > MothurOut::parseClasses(string classes){
	try {
        map<string, vector<string> > parts;
        
        //treatment<Early|Late>-age<young|old>
        vector<string> pieces; splitAtDash(classes, pieces); // -> treatment<Early|Late>, age<young|old>
        
        for (int i = 0; i < pieces.size(); i++) {
            string category = ""; string value = "";
            bool foundOpen = false;
            for (int j = 0; j < pieces[i].length(); j++) {
                if (control_pressed) { return parts; }
                
                if (pieces[i][j] == '<')        { foundOpen = true;         }
                else if (pieces[i][j] == '>')   { j += pieces[i].length();  }
                else {
                    if (!foundOpen) { category += pieces[i][j]; }
                    else { value += pieces[i][j]; }
                }
            }
            vector<string> values; splitAtChar(value, values, '|');
            parts[category] = values;
        }
        
        return parts;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "parseClasses");
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
            zapGremlins(fileHandle);
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
bool MothurOut::stringBlank(string input){
    try {
        for (int i = 0; i < input.length(); i++) {
            if (!isspace(input[i])) { return false; }
        }
        
        return true;
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
				string simpleCWD; simpleCWD="";
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
        
			fileHandle.open(completeFileName.c_str());
			if(!fileHandle) {
				//mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
				return 1;
			}else {
				//check for blank file
                zapGremlins(fileHandle);
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

		fileHandle.open(completeFileName.c_str());
		if(!fileHandle) {
			mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
			return 1;
		}
		else {
			//check for blank file
            zapGremlins(fileHandle);
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
int MothurOut::openInputFileBinary(string fileName, ifstream& fileHandle){
	try {
        
		//get full path name
		string completeFileName = getFullPathName(fileName);
        
		fileHandle.open(completeFileName.c_str(), ios::binary);
		if(!fileHandle) {
			mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
			return 1;
		}
		else {
			//check for blank file
            zapGremlins(fileHandle);
            gobble(fileHandle);
			if (fileHandle.eof()) { mothurOut("[ERROR]: " + completeFileName + " is blank. Please correct."); mothurOutEndLine();  }
			
			return 0;
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "openInputFileBinary");
		exit(1);
	}	
}
/***********************************************************************/
int MothurOut::openInputFileBinary(string fileName, ifstream& fileHandle, string noerror){
	try {
        
		//get full path name
		string completeFileName = getFullPathName(fileName);

		fileHandle.open(completeFileName.c_str(), ios::binary);
		if(!fileHandle) {
			//mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
			return 1;
		}
		else {
			//check for blank file
            zapGremlins(fileHandle);
			gobble(fileHandle);
			//if (fileHandle.eof()) { mothurOut("[ERROR]: " + completeFileName + " is blank. Please correct."); mothurOutEndLine();  }
			
			return 0;
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "openInputFileBinary - no error");
		exit(1);
	}
}
/***********************************************************************/
#ifdef USE_BOOST
int MothurOut::openInputFileBinary(string fileName, ifstream& file, boost::iostreams::filtering_istream& in){
    try {
        
        //get full path name
        string completeFileName = getFullPathName(fileName);
        
        file.open(completeFileName.c_str(), ios_base::in | ios_base::binary);
        
        if(!file) {
            mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
            return 1;
        }
        else {
            //check for blank file
            in.push(boost::iostreams::gzip_decompressor());
            in.push(file);
            if (file.eof()) { mothurOut("[ERROR]: " + completeFileName + " is blank. Please correct."); mothurOutEndLine();  }
            
            return 0;
        }
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "openInputFileGZBinary");
        exit(1);
    }
}
/***********************************************************************/
int MothurOut::openInputFileBinary(string fileName, ifstream& file, boost::iostreams::filtering_istream& in, string noerror){
    try {
        
        //get full path name
        string completeFileName = getFullPathName(fileName);
        
        file.open(completeFileName.c_str(), ios_base::in | ios_base::binary);
        
        if(!file) {
            return 1;
        }
        else {
            //check for blank file
            in.push(boost::iostreams::gzip_decompressor());
            in.push(file);
            return 0;
        }
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "openInputFileGZBinary - no error");
        exit(1);
    }
}
#endif
/***********************************************************************/
//results[0] = allGZ, results[1] = allNotGZ
vector<bool> MothurOut::allGZFiles(vector<string> & files){
    try {
        vector<bool> results;
        bool allGZ = true;
        bool allNOTGZ = true;
        
        for (int i = 0; i < files.size(); i++) {
            if (control_pressed) { break; }
            
            //ignore none and blank filenames
            if ((files[i] != "") || (files[i] != "NONE")) {
                if (isGZ(files[i])[1]) { allNOTGZ = false;  }
                else {  allGZ = false;  }
            }
        }
        
        if (!allGZ && !allNOTGZ) { //mixed bag
            mothurOut("[ERROR]: Cannot mix .gz and non compressed files. Please decompress your files and rerun.\n"); control_pressed = true; 
        }
        
        results.push_back(allGZ);
        results.push_back(allNOTGZ);
        
        return results;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "areGZFiles");
        exit(1);
    }
}

/***********************************************************************/
vector<bool> MothurOut::isGZ(string filename){
    try {
        vector<bool> results; results.resize(2, false);
        #ifdef USE_BOOST
        ifstream fileHandle;
        boost::iostreams::filtering_istream gzin;
        
        if (getExtension(filename) != ".gz") { return results; } // results[0] = false; results[1] = false;
        
        int ableToOpen = openInputFileBinary(filename, fileHandle, gzin, ""); //no error
        
        if (debug) { if (ableToOpen == 1) { mothurOut("[DEBUG]: unable to open gz file. \n"); } }
        
        if (ableToOpen == 1) { return results; } // results[0] = false; results[1] = false;
        else {  results[0] = true;  }
        
        char c;
        try
        {
            gzin >> c;
            results[1] = true;
        }
        catch ( boost::iostreams::gzip_error & e )
        {
            gzin.pop();
            fileHandle.close();
            return results;  // results[0] = true; results[1] = false;
        }
        fileHandle.close();
        #else
        mothurOut("[ERROR]: cannot test for gz format without enabling boost libraries.\n"); control_pressed = true;
        #endif
        return results; //results[0] = true; results[1] = true;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "isGZ");
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
/***********************************************************************/

int MothurOut::openOutputFileBinary(string fileName, ofstream& fileHandle){
	try {
        
		string completeFileName = getFullPathName(fileName);

		fileHandle.open(completeFileName.c_str(), ios::trunc | ios::binary);
		if(!fileHandle) {
			mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
			return 1;
		}
		else {
			return 0;
		}
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "openOutputFileBinary");
		exit(1);
	}	
    
}
/**************************************************************************************************/
int MothurOut::appendFiles(string temp, string filename) {
	try{
		ofstream output;
		ifstream input;
	
		//open output file in append mode
		openOutputFileBinaryAppend(filename, output);
		int ableToOpen = openInputFileBinary(temp, input, "no error");
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
int MothurOut::appendBinaryFiles(string temp, string filename) {
	try{
		ofstream output;
		ifstream input;
        
		//open output file in append mode
		openOutputFileBinaryAppend(filename, output);
		int ableToOpen = openInputFileBinary(temp, input, "no error");
		
		if (ableToOpen == 0) { //you opened it
            
            char buffer[4096];
            while (!input.eof()) {
                input.read(buffer, 4096);
                output.write(buffer, input.gcount());
            }
			input.close();
		}
		
		output.close();
		
		return ableToOpen;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "appendBinaryFiles");
		exit(1);
	}	
}
/**************************************************************************************************/
int MothurOut::appendSFFFiles(string temp, string filename) {
    try{
        ofstream output;
        ifstream input;
        int ableToOpen = 0;
        
        //open output file in append mode
        string fullFileName = getFullPathName(filename);
        
        output.open(fullFileName.c_str(), ios::app | ios::binary);
        if(!output) {
            mothurOut("[ERROR]: Could not open " + fullFileName); mothurOutEndLine();
            return 1;
        }else {
            //get full path name
            string completeFileName = getFullPathName(temp);
            
            input.open(completeFileName.c_str(), ios::binary);
            if(!input) {
                //mothurOut("[ERROR]: Could not open " + completeFileName); mothurOutEndLine();
                return 1;
            }else {
                char buffer[4096];
                while (!input.eof()) {
                    input.read(buffer, 4096);
                    output.write(buffer, input.gcount());
                }
                input.close();
            }
            output.close();
        }
        
        return ableToOpen;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "appendSFFFiles");
        exit(1);
    }	
}
/**************************************************************************************************/
int MothurOut::appendFilesWithoutHeaders(string temp, string filename) {
	try{
		ofstream output;
		ifstream input;
        
		//open output file in append mode
		openOutputFileAppend(filename, output);
		int ableToOpen = openInputFile(temp, input, "no error");
		//int ableToOpen = openInputFile(temp, input);
		
		int numLines = 0;
		if (ableToOpen == 0) { //you opened it
        
            string headers = getline(input); gobble(input);
            if (debug) { mothurOut("[DEBUG]: skipping headers " + headers +'\n'); }
            
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
			while (!input.eof()) {
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
            ofstream output2;
			openInputFile(tempOutfile, input2);
			openOutputFile(outfile, output2);
		
            while (!input2.eof()) {
				input2 >> dist >> firstName >> secondName;
				output2 << firstName << '\t' << secondName << '\t' << dist << endl;
				gobble(input2);
			}
			input2.close();
			output2.close();
		
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
vector<unsigned long long> MothurOut::setFilePosFasta(string filename, long long& num) {
	try {
			vector<unsigned long long> positions;
			ifstream inFASTA;
			//openInputFileBinary(filename, inFASTA);
            string completeFileName = getFullPathName(filename);
            inFASTA.open(completeFileName.c_str(), ios::binary);
						
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
					if (debug) { mothurOut("[DEBUG]: numSeqs = " + toString(positions.size()) +  " count = " + toString(count) + ".\n"); }
				}
			}
			inFASTA.close();
		
			num = positions.size();
            if (debug) { mothurOut("[DEBUG]: num = " + toString(num) + ".\n"); }
			FILE * pFile;
			unsigned long long size;
		
			//get num bytes in file
			pFile = fopen (completeFileName.c_str(),"rb");
			if (pFile==NULL) perror ("Error opening file");
			else{
				fseek (pFile, 0, SEEK_END);
				size=ftell (pFile);
				fclose (pFile);
			}
			
			/*unsigned long long size = positions[(positions.size()-1)];
			ifstream in;
			openInputFile(filename, in);
			
			in.seekg(size);
		
			while(in.get()){
				if(in.eof())		{	break;	}
				else				{	size++;	}
			}
			in.close();*/
        
            if (debug) { mothurOut("[DEBUG]: size = " + toString(size) + ".\n"); }
        
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
vector<unsigned long long> MothurOut::setFilePosFasta(string filename, long long& num, char delim) {
    try {
        vector<unsigned long long> positions;
        ifstream inFASTA;
        string completeFileName = getFullPathName(filename);
        inFASTA.open(completeFileName.c_str(), ios::binary);
        int nameLine = 2;
        if (delim == '@') { nameLine = 4; }
        else if (delim == '>') { nameLine = 2; }
        else { mothurOut("[ERROR]: unknown file deliminator, quitting.\n"); control_pressed = true; }
        
        unsigned long long count = 0;
        long long numLines = 0;
        while(!inFASTA.eof()){
            char c = inFASTA.get(); count++;
            string input = ""; input += c;
            while ((c != '\n') && (c != '\r') && (c != '\f') && (c != EOF)) {
                c = inFASTA.get(); count++;
                input += c;
            }
            numLines++;
            //gobble
            while(isspace(c=inFASTA.get()))		{ input += c; count++;}
            if(!inFASTA.eof()) { inFASTA.putback(c); count--;  }

            
            if (input.length() != 0) {
                if((input[0] == delim) && (((numLines-1)%nameLine) == 0)){ //this is a name line
                    //mothurOut(input + '\t' + toString(count+numLines-input.length()) + '\n');// << endl;
                    positions.push_back(count+numLines-input.length());
                    if (debug) { mothurOut("[DEBUG]: numSeqs = " + toString(positions.size()) +  " count = " + toString(count) + input + ".\n"); }
                }else if (int(c) == -1) { break; }
                else {
                    input = "";
                }
            }
        }
        inFASTA.close();
        
        num = positions.size();
        if (debug) { mothurOut("[DEBUG]: num = " + toString(num) + ".\n"); }
        FILE * pFile;
        unsigned long long size;
        
        //get num bytes in file
        pFile = fopen (completeFileName.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }
        
        if (debug) { mothurOut("[DEBUG]: size = " + toString(size) + ".\n"); }
        
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
vector<unsigned long long> MothurOut::setFilePosFasta(string filename, int& num) {
    try {
        vector<unsigned long long> positions;
        ifstream inFASTA;
        //openInputFile(filename, inFASTA);
        string completeFileName = getFullPathName(filename);
        inFASTA.open(completeFileName.c_str(), ios::binary);
        
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
                if (debug) { mothurOut("[DEBUG]: numSeqs = " + toString(positions.size()) +  " count = " + toString(count) + ".\n"); }
            }
        }
        inFASTA.close();
        
        num = positions.size();
        if (debug) { mothurOut("[DEBUG]: num = " + toString(num) + ".\n"); }
        FILE * pFile;
        unsigned long long size;
        
        //get num bytes in file
        pFile = fopen (completeFileName.c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }
        
        /*unsigned long long size = positions[(positions.size()-1)];
         ifstream in;
         openInputFile(filename, in);
         
         in.seekg(size);
         
         while(in.get()){
         if(in.eof())		{	break;	}
         else				{	size++;	}
         }
         in.close();*/
        
        if (debug) { mothurOut("[DEBUG]: size = " + toString(size) + ".\n"); }
        
        positions.push_back(size);
        positions[0] = 0;
        
        return positions;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "setFilePosFasta");
        exit(1);
    }
}

//**********************************************************************************************************************
vector<consTax> MothurOut::readConsTax(string inputfile){
	try {
		
        vector<consTax> taxes;
        
        ifstream in;
        openInputFile(inputfile, in);
        
        //read headers
        getline(in);
        
        while (!in.eof()) {
            
            if (control_pressed) { break; }
            
            string otu = ""; string tax = "unknown";
            int size = 0;
            
            in >> otu; gobble(in);
            in >> size; gobble(in);
            tax = getline(in); gobble(in);
            
            consTax temp(otu, tax, size);
            taxes.push_back(temp);
        }
        in.close();
        
        return taxes;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "readConsTax");
		exit(1);
	}
}
//**********************************************************************************************************************
int MothurOut::readConsTax(string inputfile, map<int, consTax2>& taxes){
	try {
        ifstream in;
        openInputFile(inputfile, in);
        
        //read headers
        getline(in);
        
        while (!in.eof()) {
            
            if (control_pressed) { break; }
            
            string otu = ""; string tax = "unknown";
            int size = 0;
            
            in >> otu; gobble(in);
            in >> size; gobble(in);
            tax = getline(in); gobble(in);
            
            consTax2 temp(otu, tax, size);
            string simpleBin = getSimpleLabel(otu);
            int bin;
            convert(simpleBin, bin);
            taxes[bin] = temp;
        }
        in.close();
        
        return 0;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "readConsTax");
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
			openInputFileBinary(filename, in);
		
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
vector<unsigned long long> MothurOut::setFilePosEachLine(string filename, unsigned long long& num) {
    try {
        filename = getFullPathName(filename);
        
        vector<unsigned long long> positions;
        ifstream in;
        //openInputFile(filename, in);
        openInputFileBinary(filename, in);
        
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
        
        if (proc > 1) {
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

vector<unsigned long long> MothurOut::divideFile(string filename, int& proc, char delimChar) {
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
        
        char secondaryDelim = '>';
        if (delimChar == '@') { secondaryDelim = '+'; }
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        
        //estimate file breaks
        unsigned long long chunkSize = 0;
        chunkSize = size / proc;
        
        //file to small to divide by processors
        if (chunkSize == 0)  {  proc = 1;	filePos.push_back(size); return filePos;	}
        
        //for each process seekg to closest file break and search for next delimChar char. make that the filebreak
        for (int i = 0; i < proc; i++) {
            unsigned long long spot = (i+1) * chunkSize;
            
            ifstream in;
            openInputFile(filename, in);
            in.seekg(spot);
            
            getline(in); //get to end of line in case you jump into middle of line where the delim char happens to fall.
            
            //look for next delimChar
            unsigned long long newSpot = spot;
            while (!in.eof()) {
                char c = in.get();
                string input = ""; input += c;
                while ((c != '\n') && (c != '\r') && (c != '\f') && (c != EOF)) {
                    c = in.get();
                    input += c;
                }
                
                if (input.length() != 0) {
                    if(input[0] == delimChar){ //this is a potential name line
                        newSpot = in.tellg();
                        newSpot -=input.length();
                        //get two lines and look for secondary delim
                        //inf a fasta file this would be a new sequence, in fastq it will be the + line, if this was a nameline.
                        getline(in); gobble(in);
                        if (!in.eof()) {
                            string secondInput = getline(in); gobble(in);
                            if (debug) { mothurOut("[DEBUG]: input= " + input + "\n secondaryInput = " + secondInput + "\n"); }
                            if (secondInput[0] == secondaryDelim) { break; } //yes, it was a nameline so stop
                            else { input = ""; gobble(in); } //nope it was a delim at the beginning of a non nameline, keep looking.
                        }
                    }else if (int(c) == -1) { break; }
                    else {  input = ""; gobble(in); }
                }
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

vector<unsigned long long> MothurOut::divideFilePerLine(string filename, int& proc) {
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
			
			//look for next line break
			unsigned long long newSpot = spot;
			while (!in.eof()) {
                char c = in.get();
				
				if ((c == '\n') || (c == '\r') || (c == '\f'))	{ gobble(in); newSpot = in.tellg(); break; }
                else if (int(c) == -1) { break; }
            }
            
			//there was not another line before the end of the file
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
/***********************************************************************/
int MothurOut::splitWhiteSpace(string input, vector<float>& pieces, int index){
    try {
        pieces.clear();
        string rest = "";
        int count = 0;
        
        for (int i = 0; i < input.length(); i++) {
            if (!isspace(input[i]))  { rest += input[i];  }
            else {
                if (rest != "") { float tdist; mothurConvert(rest, tdist); pieces.push_back(tdist); count++; rest = ""; }
                while (i < input.length()) {  //gobble white space
                    if (isspace(input[i])) { i++; }
                    else { rest = input[i];  break; } //cout << "next piece buffer = " << nextPiece << endl;
                }
                if (count > index) { return 0; }
            }
        }
        
        if (rest != "") { float tdist; mothurConvert(rest, tdist); count++; pieces.push_back(tdist); }
        
        return 0;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "splitWhiteSpace");
        exit(1);
    }
}
/***********************************************************************/
vector<string> MothurOut::splitWhiteSpaceWithQuotes(string input){
	try {
        vector<string> pieces;
        string rest = "";
        
        int pos = input.find('\'');
        int pos2 = input.find('\"');
        
        if ((pos == string::npos) && (pos2 == string::npos)) { return splitWhiteSpace(input); } //no quotes to worry about
        else {
            for (int i = 0; i < input.length(); i++) {
                if ((input[i] == '\'') || (input[i] == '\"') || (rest == "\'") || (rest == "\"")) { //grab everything til end or next ' or "
                    rest += input[i];
                    for (int j = i+1; j < input.length(); j++) {
                        if ((input[j] == '\'') || (input[j] == '\"')) {  //then quit
                            rest += input[j];
                            i = j+1;
                            j+=input.length();
                        }else { rest += input[j]; }
                    }
                }else if (!isspace(input[i]))  { rest += input[i];  }
                else {
                    if (rest != "") { pieces.push_back(rest);  rest = ""; }
                    while (i < input.length()) {  //gobble white space
                        if (isspace(input[i])) { i++; }
                        else { rest = input[i];  break; } //cout << "next piece buffer = " << nextPiece << endl;
                    } 
                }
            }
            
            if (rest != "") { pieces.push_back(rest); }
        }
        return pieces;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitWhiteSpace");
		exit(1);
	}
}
//**********************************************************************************************************************
int MothurOut::readTax(string taxfile, map<string, string>& taxMap, bool removeConfidence) {
	try {
        //open input file
		ifstream in;
		openInputFile(taxfile, in);
        
        bool error = false;
        string name, taxonomy;
        
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in >> name; gobble(in);
            taxonomy = getline(in); gobble(in);
            
            checkName(name);
            //are there confidence scores, if so remove them
            if (removeConfidence) {  if (taxonomy.find_first_of('(') != -1) {  removeConfidences(taxonomy);	} }
            map<string, string>::iterator itTax = taxMap.find(name);
            
            if(itTax == taxMap.end()) {
                bool ignore = false;
                if (taxonomy != "") { if (taxonomy[taxonomy.length()-1] != ';') { mothurOut("[ERROR]: " + name + " is missing the final ';', ignoring.\n"); ignore=true; }
                }
                if (!ignore) { taxMap[name] = taxonomy; }
                if (debug) {  mothurOut("[DEBUG]: name = '" + name + "' tax = '" + taxonomy + "'\n");  }
            }else {
                mothurOut("[ERROR]: " + name + " is already in your taxonomy file, names must be unique./n"); error = true;
            }

        }
		in.close();
        
        if (error) { control_pressed = true; }
        if (debug) {  mothurOut("[DEBUG]: numSeqs saved = '" + toString(taxMap.size()) + "'\n"); }
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
                    checkName(firstCol);
                    checkName(secondCol);
                    
                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = firstCol;  }
                    pairDone = false; 
                }
            }
		}
		in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    
                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {   nameMap[theseNames[i]] = firstCol;  }
                    pairDone = false; 
                }
            }  
        }
		
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
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[secondCol] = firstCol;
                    pairDone = false; 
                }
            }
		}
		in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[secondCol] = firstCol;
                    pairDone = false; 
                }
            } 
        }
		
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
                    checkName(firstCol);
                    checkName(secondCol);
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
		
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    //parse names into vector
                    vector<string> theseNames;
                    splitAtComma(secondCol, theseNames);
                    for (int i = 0; i < theseNames.size(); i++) {  nameMap[theseNames[i]] = firstCol;  }
                    nameCount[firstCol] = theseNames.size();
                    pairDone = false; 
                }
            }

        }
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
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[firstCol] = secondCol; pairDone = false; }
            }
		}
		in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    nameMap[firstCol] = secondCol; pairDone = false; }
            }
        }
		
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
                    checkName(firstCol);
                    checkName(secondCol);
                    vector<string> temp;
                    splitAtComma(secondCol, temp);
                    nameMap[firstCol] = temp;
                    pairDone = false;  
                } 
            }
		}
		in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    vector<string> temp;
                    splitAtComma(secondCol, temp);
                    nameMap[firstCol] = temp;
                    pairDone = false;  
                } 
            }
        }
        
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
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    if (debug) { mothurOut("[DEBUG]: " + firstCol + ", " + toString(num) + ".\n"); }
                    pairDone = false;  
                } 
            }
		}
        in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;  
                } 
            }
        }
		
		return nameMap;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
map<string, int> MothurOut::readNames(string namefile, unsigned long int& numSeqs) { 
	try {
		map<string, int> nameMap;
        numSeqs = 0;
		
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
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;  
                    numSeqs += num;
                } 
            }
		}
        in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
                    int num = getNumNames(secondCol);
                    nameMap[firstCol] = num;
                    pairDone = false;  
                    numSeqs += num;
                } 
            }
        }
		
		return nameMap;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readNames");
		exit(1);
	}
}
/************************************************************/
int MothurOut::checkName(string& name) {
    try {
        if (modifyNames) {
            for (int i = 0; i < name.length(); i++) {
                if (name[i] == ':') { name[i] = '_'; changedSeqNames = true; }
            }
        }
        return 0;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "checkName");
		exit(1);
	}
}
/************************************************************/
bool MothurOut::checkGroupName(string name) {
    try {
        
        bool goodName = true;
        for (int i = 0; i < name.length(); i++) {
            if (name[i] == ':') {  goodName = false; break;  }
            else if (name[i] == '-') {  goodName = false; break;  }
            else if (name[i] == '/') {  goodName = false; break;  }
        }
        
        if (!goodName) {
            mothurOut("\n[WARNING]: group " + name + " contains illegal characters in the name. Group names should not include :, -, or / characters.  The ':' character is a special character used in trees. Using ':' will result in your tree being unreadable by tree reading software.  The '-' character is a special character used by mothur to parse group names.  Using the '-' character will prevent you from selecting groups. The '/' character will created unreadable filenames when mothur includes the group in an output filename.  \n\n");
        }
        
        return goodName;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "checkGroupName");
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
                    checkName(firstCol);
                    checkName(secondCol);
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
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  firstCol = pieces[i]; columnOne=false; }
                else  { secondCol = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    checkName(firstCol);
                    checkName(secondCol);
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
		int ableToOpen = openInputFile(accnosfile, in, "");
        if (ableToOpen == 1) {  mothurOut("[ERROR]: Could not open " + accnosfile); mothurOutEndLine(); return names; }
		string name;
		
        string rest = "";
        char buffer[4096];
        
        unsigned long long count = 0;
		while (!in.eof()) {
			if (control_pressed) { break; }
			
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]);
                names.insert(pieces[i]);
                count++;
            }
        }
		in.close();	
		
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.insert(pieces[i]);  count++; }
        }
        
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
            
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }
		in.close();	
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "readAccnos");
		exit(1);
	}
}
//**********************************************************************************************************************
int MothurOut::readAccnos(string accnosfile, vector<string>& names, string noerror){
    try {
        names.clear();
        ifstream in;
        openInputFile(accnosfile, in, noerror);
        string name;
        
        string rest = "";
        char buffer[4096];
        
        while (!in.eof()) {
            if (control_pressed) { break; }
            
            in.read(buffer, 4096);
            vector<string> pieces = splitWhiteSpace(rest, buffer, in.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }
        in.close();
        
        if (rest != "") {
            vector<string> pieces = splitWhiteSpace(rest);
            for (int i = 0; i < pieces.size(); i++) {  checkName(pieces[i]); names.push_back(pieces[i]);  }
        }
        
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
string MothurOut::getSimpleLabel(string label){
	try {
		string simple = "";
        
        //remove OTU or phylo tag
        string newLabel1 = "";
        for (int i = 0; i < label.length(); i++) {
            if(label[i]>47 && label[i]<58) { //is a digit
                newLabel1 += label[i];
            }
        }
        
        int num1;
        
        if (debug) { mothurOut("[DEBUG]: " + newLabel1 + "\n"); }
            
        mothurConvert(newLabel1, num1);
        
        simple = toString(num1);
        
		return simple;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "getSimpleLabel");
		exit(1);
	}
}
/***********************************************************************/
string MothurOut::mothurGetpid(int threadID){
	try {
        
        string pid = "";
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        
		pid += toString(getpid()); if(debug) { mothurOut("[DEBUG]: " + pid + "\n"); }
        //remove any weird chars
        string pid1 = "";
        for (int i = 0; i < pid.length(); i++) {
            if(pid[i]>47 && pid[i]<58) { //is a digit
                pid1 += pid[i];
            }
        }
        pid = pid1;
#else
		pid += toString(threadID);
#endif
		return pid;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "mothurGetpid");
		exit(1);
	}
}

/***********************************************************************/

bool MothurOut::isLabelEquivalent(string label1,  string label2){
	try {
		bool same = false;
        
        //remove OTU or phylo tag
        string newLabel1 = "";
        for (int i = 0; i < label1.length(); i++) {
            if(label1[i]>47 && label1[i]<58) { //is a digit
                newLabel1 += label1[i];
            }
        }
        
        string newLabel2 = "";
        for (int i = 0; i < label2.length(); i++) {
            if(label2[i]>47 && label2[i]<58) { //is a digit
                newLabel2 += label2[i];
            }
        }
        
        int num1, num2;
        mothurConvert(newLabel1, num1);
        mothurConvert(newLabel2, num2);
        
        if (num1 == num2) { same = true; }
		
		return same;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "isLabelEquivalent");
		exit(1);
	}
}
//**********************************************************************************************************************
bool MothurOut::isSubset(vector<string> bigset, vector<string> subset) {
	try {
		
        
		if (subset.size() > bigset.size()) { return false;  }
		
		//check if each guy in subset is also in bigset
		for (int i = 0; i < subset.size(); i++) {
			bool match = false;
			for (int j = 0; j < bigset.size(); j++) {
				if (subset[i] == bigset[j]) { match = true; break; }
			}
			
			//you have a guy in subset that had no match in bigset
			if (match == false) { return false; }
		}
		
		return true;
        
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "isSubset");
		exit(1);
	}
}
/***********************************************************************/
int MothurOut::mothurRemove(string filename){
	try {
		filename = getFullPathName(filename);
		int error = remove(filename.c_str());
        if (debug) {
            mothurOut("[DEBUG]: removed " + filename + "\n");
        }
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
		errorOut(e, "MothurOut", "mothurConvert-int");
		exit(1);
	}
}
/***********************************************************************/
bool MothurOut::mothurConvert(char item, int& num){
    try {
        bool error = false;
        
        if (isdigit(item)) {
            string mystring; mothurConvert(item, mystring);
            mothurConvert(mystring, num);
        }else {
            num = 0;
            error = true;
            mothurOut("[ERROR]: cannot convert " + toString(item) + " to an integer."); mothurOutEndLine();
            commandInputsConvertError = true;
        }
        
        return error;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "mothurConvert-int");
        exit(1);
    }
}
/***********************************************************************/
bool MothurOut::mothurConvert(char item, string& output){
    try {
        
        stringstream ss;
        ss << item;
        ss >> output;
        return true;
        
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "mothurConvert-char");
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
		errorOut(e, "MothurOut", "mothurConvert-intDist");
		exit(1);
	}
}
/**************************************************************************************************/
string MothurOut::addUnclassifieds(string tax, int maxlevel, bool probs) {
    try{
        string newTax, taxon;
        
        string savedTax = tax;
        vector<string> taxons; splitAtChar(tax, taxons, ';'); taxons.pop_back();
        vector<int> confidences;
        
        if (taxons.size() == maxlevel) { return savedTax; }
        
        int index = 0;
        int confidence = 0;
        int level = 1;
        for (int i = 0; i < taxons.size(); i++) {
            index = i;
            string thisTax = taxons[i]+";";
            confidence = removeConfidences(thisTax);
            confidences.push_back(confidence);
    
            if (thisTax == "unclassified;"){ index--; break; }
            else{ newTax += taxons[i] + ";";  }
        }
        level = index+1;
        
        string thisTax = taxons[index]+";";
        
        removeConfidences(thisTax);
        taxon = thisTax.substr(0, thisTax.length()-1);
        
        string cTax = "";
        if (probs)  { cTax = taxon + "_unclassified(" + toString(confidences[index]) + ");";     }
        else        { cTax = taxon + "_unclassified;";          }
        
        //add "unclassified" until you reach maxLevel
        while (level < maxlevel) {
            newTax += cTax;
            level++;
        }
        
        return newTax;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "addUnclassifieds");
        exit(1);
    }
}

/***********************************************************************/
bool MothurOut::isNumeric1(string stringToCheck){
	try {
		bool numeric = false;
		
        if (stringToCheck == "") { numeric = false;  }
        else if(stringToCheck.find_first_not_of("0123456789.-") == string::npos) { numeric = true; }
			
		return numeric;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "isNumeric1");
		exit(1);
	}
	
}
/***********************************************************************/
bool MothurOut::allSpaces(string stringToCheck){
    try {
        
        for (int i = 0; i < stringToCheck.length(); i++) {
            char c = stringToCheck[i];
            if (!isspace(c)) { return false; }
        }
        
        return true;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "isNumeric1");
        exit(1);
    }
    
}
/***********************************************************************/
bool MothurOut::isInteger(string stringToCheck){
    try {
        bool isInt = false;
        
        if(stringToCheck.find_first_not_of("0123456789-") == string::npos) { isInt = true; }
        
        return isInt;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "isInteger");
        exit(1);
    }
    
}
/***********************************************************************/
bool MothurOut::containsAlphas(string stringToCheck){
    try {
        bool containsAlpha = false;
        
        if(stringToCheck.find_first_of("AaBbCcDdEeFfGgHhIiJjKkLlMmNnOopPQqRrSsTtUuVvWwXxYyZz") != string::npos) { containsAlpha = true; }
        
        return containsAlpha;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "containsAlphas");
        exit(1);
    }
    
}
/***********************************************************************/
bool MothurOut::isAllAlphas(string stringToCheck){
    try {
        bool allAlphas = true;
        
        if(stringToCheck.find_first_not_of("AaBbCcDdEeFfGgHhIiJjKkLlMmNnOopPQqRrSsTtUuVvWwXxYyZz") != string::npos) { allAlphas = false; }
        
        return allAlphas;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "isAllAlphas");
        exit(1);
    }
    
}
/***********************************************************************/
bool MothurOut::isAllAlphaNumerics(string stringToCheck){
    try {
        bool allAlphaNumerics = true;
        
        if(stringToCheck.find_first_not_of("AaBbCcDdEeFfGgHhIiJjKkLlMmNnOopPQqRrSsTtUuVvWwXxYyZz0123456789") != string::npos) { allAlphaNumerics = false; }
        
        return allAlphaNumerics;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "isAllAlphas");
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
		errorOut(e, "MothurOut", "mothurConvert-float");
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
		errorOut(e, "MothurOut", "mothurConvert-double");
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
string  MothurOut::findEdianness() {
    try {
        // find real endian type
        string endianType = "unknown";
        int num = 1;
        if(*(char *)&num == 1)
        {
            endianType = "LITTLE_ENDIAN";
        }
        else
        {
            endianType = "BIG_ENDIAN";
        }
        return endianType;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "findEdianness");
		exit(1);
	}
}
/***********************************************************************/
double  MothurOut::median(vector<double> x) {
    try {
        double value = 0.0;
        
        if (x.size() == 0) { } //error
        else {
            //For example, if a < b < c, then the median of the list {a, b, c} is b, and, if a < b < c < d, then the median of the list {a, b, c, d} is the mean of b and c; i.e., it is (b + c)/2.
            sort(x.begin(), x.end());
            //is x.size even?
            if ((x.size()%2) == 0) { //size() is even. median = average of 2 midpoints
                int midIndex1 = (x.size()/2)-1;
                int midIndex2 = (x.size()/2);
                value = (x[midIndex1]+ x[midIndex2]) / 2.0;
            }else { 
                int midIndex = (x.size()/2);
                value = x[midIndex];
            }
        }
        return value;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "median");
		exit(1);
	}
}
/***********************************************************************/
int  MothurOut::median(vector<int> x) {
    try {
        double value = 0;
        
        if (x.size() == 0) { } //error
        else {
            //For example, if a < b < c, then the median of the list {a, b, c} is b, and, if a < b < c < d, then the median of the list {a, b, c, d} is the mean of b and c; i.e., it is (b + c)/2.
            sort(x.begin(), x.end());
            //is x.size even?
            if ((x.size()%2) == 0) { //size() is even. median = average of 2 midpoints
                int midIndex1 = (x.size()/2)-1;
                int midIndex2 = (x.size()/2);
                value = (x[midIndex1]+ x[midIndex2]) / 2.0;
            }else {
                int midIndex = (x.size()/2);
                value = x[midIndex];
            }
        }
        return (int) value;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "median - int");
        exit(1);
    }
}
/***********************************************************************/
int  MothurOut::average(vector<int> x) {
    try {
        int value = 0;
        
        for (int i = 0; i < x.size(); i++) {
            if (control_pressed) { break; }
            value += x[i];
        }
        
        return ((int) value / x.size());
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "average - int");
        exit(1);
    }
}
/***********************************************************************/
int  MothurOut::sum(vector<int> x) {
    try {
        int value = 0;
        
        for (int i = 0; i < x.size(); i++) {
            if (control_pressed) { break; }
            value += x[i];
        }
        
        return value;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "sum - int");
        exit(1);
    }
}
/***********************************************************************/
double  MothurOut::sum(vector<double> x) {
    try {
        int value = 0;
        
        for (int i = 0; i < x.size(); i++) {
            if (control_pressed) { break; }
            value += x[i];
        }
        
        return value;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "sum - double");
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
bool MothurOut::checkLocations(string& filename, string inputDir){
	try {
		filename = getFullPathName(filename);
        
        int ableToOpen;
        ifstream in;
        ableToOpen = openInputFile(filename, in, "noerror");
        in.close();
        
        //if you can't open it, try input location
        if (ableToOpen == 1) {
            if (inputDir != "") { //default path is set
                string tryPath = inputDir + getSimpleName(filename);
                mothurOut("Unable to open " + filename + ". Trying input directory " + tryPath); mothurOutEndLine();
                ifstream in2;
                ableToOpen = openInputFile(tryPath, in2, "noerror");
                in2.close();
                filename = tryPath;
            }
        }
        
        //if you can't open it, try default location
        if (ableToOpen == 1) {
            if (getDefaultPath() != "") { //default path is set
                string tryPath = getDefaultPath() + getSimpleName(filename);
                mothurOut("Unable to open " + filename + ". Trying default " + tryPath); mothurOutEndLine();
                ifstream in2;
                ableToOpen = openInputFile(tryPath, in2, "noerror");
                in2.close();
                filename = tryPath;
            }
        }
        
        //if you can't open it its not in current working directory or inputDir, try mothur excutable location
        if (ableToOpen == 1) {
            string exepath = mothurProgramPath;
            string tempPath = exepath;
            for (int i = 0; i < exepath.length(); i++) { tempPath[i] = tolower(exepath[i]); }
            exepath = exepath.substr(0, (tempPath.find_last_of('m')));
            
            string tryPath = getFullPathName(exepath) + getSimpleName(filename);
            mothurOut("Unable to open " + filename + ". Trying mothur's executable location " + tryPath); mothurOutEndLine();
            ifstream in2;
            ableToOpen = openInputFile(tryPath, in2, "noerror");
            in2.close();
            filename = tryPath;
        }
        
        if (ableToOpen == 1) { mothurOut("Unable to open " + filename + "."); mothurOutEndLine(); return false;  }
        
        return true;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "checkLocations");
		exit(1);
	}
}
/***********************************************************************/

//This function parses the estimator options and puts them in a vector
void MothurOut::splitAtChar(string& estim, vector<string>& container, char symbol) {
	try {
        
        if (symbol == '-') { splitAtDash(estim, container); return; }
        
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
		bool prevEscape = false;
        
        for(int i=0;i<estimLength;i++){
            if(estim[i] == '-'){
                if (prevEscape) {  individual += estim[i]; prevEscape = false;  } //add in dash because it was escaped.
                else {
                    container.push_back(individual);
                    individual = "";
                }
            }else if(estim[i] == '\\'){
                if (i < estimLength-1) { 
                    if (estim[i+1] == '-') { prevEscape=true; }  //are you a backslash before a dash, if yes ignore
                    else { individual += estim[i]; prevEscape = false;  } //if no, add in
                }else { individual += estim[i]; }
            }else {
                individual += estim[i];
            }
        }
        

        
		container.push_back(individual);
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
		bool prevEscape = false;
        
        for(int i=0;i<estimLength;i++){
            if(estim[i] == '-'){
                if (prevEscape) {  individual += estim[i]; prevEscape = false;  } //add in dash because it was escaped.
                else {
                    container.insert(individual);
                    individual = "";
                }
            }else if(estim[i] == '\\'){
                if (i < estimLength-1) { 
                    if (estim[i+1] == '-') { prevEscape=true; }  //are you a backslash before a dash, if yes ignore
                    else { individual += estim[i]; prevEscape = false;  } //if no, add in
                }else { individual += estim[i]; }
            }else {
                individual += estim[i];
            }
        }
        container.insert(individual);
        
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
		string individual = "";
		int lineNum;
		int estimLength = estim.size();
		bool prevEscape = false;
        
        for(int i=0;i<estimLength;i++){
            if(estim[i] == '-'){
                if (prevEscape) {  individual += estim[i]; prevEscape = false;  } //add in dash because it was escaped.
                else {
                    convert(individual, lineNum); //convert the string to int
                    container.insert(lineNum);
                    individual = "";
                }
            }else if(estim[i] == '\\'){
                if (i < estimLength-1) { 
                    if (estim[i+1] == '-') { prevEscape=true; }  //are you a backslash before a dash, if yes ignore
                    else { individual += estim[i]; prevEscape = false;  } //if no, add in
                }else { individual += estim[i]; }
            }else {
                individual += estim[i];
            }
        }
        
		convert(individual, lineNum); //convert the string to int
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
        
        string individual = "";
        int estimLength = prefix.size();
        for(int i=0;i<estimLength;i++){
            if(prefix[i] == c){
                suffix = prefix.substr(i+1);
                prefix = individual;
                break;
            }
            else{
                individual += prefix[i];
            }
        }
        
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "splitAtChar");
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
		}else {  suffix = "";  }

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

bool MothurOut::inUsersGroups(vector<int> set, vector< vector<int> > sets) {
	try {
		for (int i = 0; i < sets.size(); i++) {
			if (set == sets[i]) { return true; }
		}
		return false;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "inUsersGroups");
		exit(1);
	}	
}
/**************************************************************************************************/

bool MothurOut::inUsersGroups(int groupname, vector<int> Groups) {
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
/**************************************************************************************************/
//removes entries that are only white space
int MothurOut::removeBlanks(vector<string>& tempVector) {
	try {
		vector<string> newVector;
		for (int i = 0; i < tempVector.size(); i++) {
            bool isBlank = true;
            for (int j = 0; j < tempVector[i].length(); j++) {
                if (!isspace(tempVector[i][j])) { isBlank = false; j+= tempVector[i].length(); } //contains non space chars, break out and save
            }
            if (!isBlank) { newVector.push_back(tempVector[i]); }
        }
        tempVector = newVector;
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "removeBlanks");
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
int MothurOut::getTimeStamp(string filename) {
    try {
        int timeStamp = 0;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        struct stat st;
        int errorCode = stat (filename.c_str(), &st);
        if (errorCode != 0) {
            mothurOut("[ERROR]: Can't find timestamp for " + filename + "\n"); control_pressed = true;
        }else {
            timeStamp = st.st_mtime;
        }
#else
        HANDLE hFile;
        
        hFile = CreateFile(filename.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL,
                           OPEN_EXISTING, 0, NULL);
        
        if(hFile == INVALID_HANDLE_VALUE) {
            mothurOut("[ERROR]: Can't find timestamp for " + filename + "\n"); control_pressed = true;
            return timeStamp;
        }
            
        FILETIME ftCreate, ftAccess, ftWrite;
        SYSTEMTIME stUTC;
        DWORD dwRet;
        
        // Retrieve the file times for the file.
        bool success = GetFileTime(hFile, &ftCreate, &ftAccess, &ftWrite);
        
        if (success) {
            FileTimeToSystemTime(&ftWrite, &stUTC);
            
            tm time;
            time.tm_sec = stUTC.wSecond;
            time.tm_min = stUTC.wMinute;
            time.tm_hour = stUTC.wHour;
            time.tm_mday = stUTC.wDay;
            time.tm_mon = stUTC.wMonth - 1;
            time.tm_year = stUTC.wYear - 1900;
            time.tm_isdst = -1;
            time_t t = mktime(&time);
            
            timeStamp = t;
        }
        else { mothurOut("[ERROR]: Can't find timestamp for " + filename + "\n"); control_pressed = true; }
        
#endif
    
        return timeStamp;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "getTimeStamp");
        exit(1);
    }
}
/**************************************************************************************************/
vector<double> MothurOut::getAverages(vector< vector<double> >& dists) {
	try{
        vector<double> averages; //averages.resize(numComp, 0.0);
        for (int i = 0; i < dists[0].size(); i++) { averages.push_back(0.0); }
      
        for (int thisIter = 0; thisIter < dists.size(); thisIter++) {
            for (int i = 0; i < dists[thisIter].size(); i++) {  
                averages[i] += dists[thisIter][i];
            }
        }
        
        //finds average.
        for (int i = 0; i < averages.size(); i++) {  averages[i] /= (double) dists.size(); }
        
        return averages;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getAverages");		
		exit(1);
	}
}
/**************************************************************************************************/
double MothurOut::getAverage(vector<double> dists) {
	try{
        double average = 0;
        
        for (int i = 0; i < dists.size(); i++) {
            average += dists[i];
        }
       
        //finds average.
        average /= (double) dists.size(); 
        
        return average;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getAverage");
		exit(1);
	}
}

/**************************************************************************************************/
vector<double> MothurOut::getStandardDeviation(vector< vector<double> >& dists) {
	try{
        
        vector<double> averages = getAverages(dists);
        
        //find standard deviation
        vector<double> stdDev; //stdDev.resize(numComp, 0.0);
        for (int i = 0; i < dists[0].size(); i++) { stdDev.push_back(0.0); }
        
        for (int thisIter = 0; thisIter < dists.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int j = 0; j < dists[thisIter].size(); j++) {
                stdDev[j] += ((dists[thisIter][j] - averages[j]) * (dists[thisIter][j] - averages[j]));
            }
        }
        for (int i = 0; i < stdDev.size(); i++) {  
            stdDev[i] /= (double) dists.size(); 
            stdDev[i] = sqrt(stdDev[i]);
        }
        
        return stdDev;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getAverages");		
		exit(1);
	}
}
/**************************************************************************************************/
vector<double> MothurOut::getStandardDeviation(vector< vector<double> >& dists, vector<double>& averages) {
	try{
        //find standard deviation
        vector<double> stdDev; //stdDev.resize(numComp, 0.0);
        for (int i = 0; i < dists[0].size(); i++) { stdDev.push_back(0.0); }
        
        for (int thisIter = 0; thisIter < dists.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int j = 0; j < dists[thisIter].size(); j++) {
                stdDev[j] += ((dists[thisIter][j] - averages[j]) * (dists[thisIter][j] - averages[j]));
            }
        }
        for (int i = 0; i < stdDev.size(); i++) {  
            stdDev[i] /= (double) dists.size(); 
            stdDev[i] = sqrt(stdDev[i]);
        }
        
        return stdDev;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getStandardDeviation");		
		exit(1);
	}
}
/**************************************************************************************************/
vector< vector<seqDist> > MothurOut::getAverages(vector< vector< vector<seqDist> > >& calcDistsTotals, string mode) {
	try{
        
        vector< vector<seqDist>  > calcAverages; //calcAverages.resize(calcDistsTotals[0].size()); 
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            //calcAverages[i].resize(calcDistsTotals[0][i].size());
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            calcAverages.push_back(temp);
        }
        
        if (mode == "average") {
            for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //sum all groups dists for each calculator
                for (int i = 0; i < calcAverages.size(); i++) {  //initialize sums to zero.
                    for (int j = 0; j < calcAverages[i].size(); j++) {
                        calcAverages[i][j].dist += calcDistsTotals[thisIter][i][j].dist;
                    }
                }
            }
            
            for (int i = 0; i < calcAverages.size(); i++) {  //finds average.
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    calcAverages[i][j].dist /= (float) calcDistsTotals.size();
                }
            }
        }else { //find median
            for (int i = 0; i < calcAverages.size(); i++) { //for each calc
                for (int j = 0; j < calcAverages[i].size(); j++) {  //for each comparison
                    vector<double> dists;
                    for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //for each subsample
                        dists.push_back(calcDistsTotals[thisIter][i][j].dist);
                    }
                    sort(dists.begin(), dists.end());
                    calcAverages[i][j].dist = dists[(calcDistsTotals.size()/2)];
                }
            }
        }

        return calcAverages;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getAverages");		
		exit(1);
	}
}
/**************************************************************************************************/
vector< vector<seqDist> > MothurOut::getAverages(vector< vector< vector<seqDist> > >& calcDistsTotals) {
	try{
        
        vector< vector<seqDist>  > calcAverages; //calcAverages.resize(calcDistsTotals[0].size()); 
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            //calcAverages[i].resize(calcDistsTotals[0][i].size());
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            calcAverages.push_back(temp);
        }
        
        
        for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //sum all groups dists for each calculator
                for (int i = 0; i < calcAverages.size(); i++) {  //initialize sums to zero.
                    for (int j = 0; j < calcAverages[i].size(); j++) {
                        calcAverages[i][j].dist += calcDistsTotals[thisIter][i][j].dist;
                    }
                }
        }
            
        for (int i = 0; i < calcAverages.size(); i++) {  //finds average.
                for (int j = 0; j < calcAverages[i].size(); j++) {
                    calcAverages[i][j].dist /= (float) calcDistsTotals.size();
                }
        }
        
        return calcAverages;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getAverages");		
		exit(1);
	}
}
/**************************************************************************************************/
vector< vector<seqDist> > MothurOut::getStandardDeviation(vector< vector< vector<seqDist> > >& calcDistsTotals) {
	try{
        
        vector< vector<seqDist> > calcAverages = getAverages(calcDistsTotals);
        
        //find standard deviation
        vector< vector<seqDist>  > stdDev;  
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            stdDev.push_back(temp);
        }
        
        for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int i = 0; i < stdDev.size(); i++) {  
                for (int j = 0; j < stdDev[i].size(); j++) {
                    stdDev[i][j].dist += ((calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist) * (calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist));
                }
            }
        }
        
        for (int i = 0; i < stdDev.size(); i++) {  //finds average.
            for (int j = 0; j < stdDev[i].size(); j++) {
                stdDev[i][j].dist /= (float) calcDistsTotals.size();
                stdDev[i][j].dist = sqrt(stdDev[i][j].dist);
            }
        }

        return stdDev;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getAverages");		
		exit(1);
	}
}
/**************************************************************************************************/
vector< vector<seqDist> > MothurOut::getStandardDeviation(vector< vector< vector<seqDist> > >& calcDistsTotals, vector< vector<seqDist> >& calcAverages) {
	try{
        //find standard deviation
        vector< vector<seqDist>  > stdDev;  
        for (int i = 0; i < calcDistsTotals[0].size(); i++) {  //initialize sums to zero.
            vector<seqDist> temp;
            for (int j = 0; j < calcDistsTotals[0][i].size(); j++) {
                seqDist tempDist;
                tempDist.seq1 = calcDistsTotals[0][i][j].seq1;
                tempDist.seq2 = calcDistsTotals[0][i][j].seq2;
                tempDist.dist = 0.0;
                temp.push_back(tempDist);
            }
            stdDev.push_back(temp);
        }
        
        for (int thisIter = 0; thisIter < calcDistsTotals.size(); thisIter++) { //compute the difference of each dist from the mean, and square the result of each
            for (int i = 0; i < stdDev.size(); i++) {  
                for (int j = 0; j < stdDev[i].size(); j++) {
                    stdDev[i][j].dist += ((calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist) * (calcDistsTotals[thisIter][i][j].dist - calcAverages[i][j].dist));
                }
            }
        }
        
        for (int i = 0; i < stdDev.size(); i++) {  //finds average.
            for (int j = 0; j < stdDev[i].size(); j++) {
                stdDev[i][j].dist /= (float) calcDistsTotals.size();
                stdDev[i][j].dist = sqrt(stdDev[i][j].dist);
            }
        }
        
        return stdDev;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getAverages");		
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
float MothurOut::removeConfidences(string& tax) {
	try {
		
		string taxon;
		string newTax = "";
        string confidenceScore = "0";
		
		while (tax.find_first_of(';') != -1) {
			
			if (control_pressed) { return 0; }
			
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'));
	
			int pos = taxon.find_last_of('(');
			if (pos != -1) {
				//is it a number?
				int pos2 = taxon.find_last_of(')');
				if (pos2 != -1) {
					string temp = taxon.substr(pos+1, (pos2-(pos+1)));
					if (isNumeric1(temp)) {
						taxon = taxon.substr(0, pos); //rip off confidence
                        confidenceScore = temp;
					}
				}
			}
			taxon += ";";
			
			tax = tax.substr(tax.find_first_of(';')+1, tax.length());
			newTax += taxon;
		}
		
		tax = newTax;
        
        float confidence = 0; mothurConvert(confidenceScore, confidence);
		
		return confidence;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "removeConfidences");
		exit(1);
	}
}
/**************************************************************************************************/
string MothurOut::removeQuotes(string tax) {
	try {
		
		string taxon;
		string newTax = "";
		
		for (int i = 0; i < tax.length(); i++) {
			
			if (control_pressed) { return newTax; }
            
            if ((tax[i] != '\'') && (tax[i] != '\"')) { newTax += tax[i]; }
			
        }
		
		return newTax;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "removeQuotes");
		exit(1);
	}
}
/**************************************************************************************************/
// function for calculating standard deviation
double MothurOut::getStandardDeviation(vector<int>& featureVector){
    try {
        //finds sum
        double average = 0; 
        for (int i = 0; i < featureVector.size(); i++) { average += featureVector[i]; }
        average /= (double) featureVector.size();
        
        //find standard deviation
        double stdDev = 0;
        for (int i = 0; i < featureVector.size(); i++) { //compute the difference of each dist from the mean, and square the result of each
            stdDev += ((featureVector[i] - average) * (featureVector[i] - average));
        }
          
        stdDev /= (double) featureVector.size(); 
        stdDev = sqrt(stdDev);
        
        return stdDev;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "getStandardDeviation");
		exit(1);
	}
}
/**************************************************************************************************/
// returns largest value in vector
double MothurOut::max(vector<double>& featureVector){
    try {
        if (featureVector.size() == 0) { mothurOut("[ERROR]: vector size = 0!\n"); control_pressed=true; return 0.0; }
        
        //finds largest
        double largest = featureVector[0];
        for (int i = 1; i < featureVector.size(); i++) {
            if (featureVector[i] > largest) { largest = featureVector[i]; }
        }
                
        return largest;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "max");
		exit(1);
	}
}
/**************************************************************************************************/
// returns smallest value in vector
double MothurOut::min(vector<double>& featureVector){
    try {
        if (featureVector.size() == 0) { mothurOut("[ERROR]: vector size = 0!\n"); control_pressed=true; return 0.0; }
        
        //finds smallest
        double smallest = featureVector[0];
        for (int i = 1; i < featureVector.size(); i++) {
            if (featureVector[i] < smallest) { smallest = featureVector[i]; }
        }
        
        return smallest;
    }
	catch(exception& e) {
		errorOut(e, "MothurOut", "min");
		exit(1);
	}
}
/**************************************************************************************************/
int MothurOut::max(int A, int B){
    try {
        
        //finds largest
        int largest = A;
        if (B > A) { largest = B; }
        
        return largest;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "max");
        exit(1);
    }
}
/**************************************************************************************************/
int MothurOut::min(int A, int B){
    try {
        
        //finds smallest
        int smallest = A;
        if (B < A) { smallest = B; }
        
        return smallest;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "min");
        exit(1);
    }
}
//**********************************************************************************************************************
int MothurOut::printVsearchFile(vector<seqPriorityNode>& nameMapCount, string filename, string tag, string tag2){
    try {
        
        sort(nameMapCount.begin(), nameMapCount.end(), compareSeqPriorityNodes);
        
        ofstream out;
        openOutputFile(filename, out);
        
        //print new file in order of
        for (int i = 0; i < nameMapCount.size(); i++) {
            if (control_pressed) {break;}
            out << ">" << nameMapCount[i].name  << tag << nameMapCount[i].numIdentical << tag2 << endl << nameMapCount[i].seq << endl;
        }
        out.close();
        
        return 0;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "printVsearchFile");
        exit(1);
    }
}
//**********************************************************************************************************************
string MothurOut::getStringFromVector(vector<string>& list, string delim){
    try {
        string result = "";
        
        if (list.size() == 0) { return result; }
        
        result = list[0];
        
        for (int i = 1; i < list.size(); i++) {
            if (control_pressed) { break;  }
            result += delim + list[i];
        }
        
        return result;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "getStringFromVector");
        exit(1);
    }
}
//**********************************************************************************************************************
string MothurOut::getStringFromVector(vector<int>& list, string delim){
    try {
        string result = "";
        
        if (list.size() == 0) { return result; }
        
        result = list[0];
        
        for (int i = 1; i < list.size(); i++) {
            if (control_pressed) { break;  }
            string temp = toString(list[i]);
            result += delim + temp;
        }
        
        return result;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "getStringFromVector");
        exit(1);
    }
}
//**********************************************************************************************************************
string MothurOut::getStringFromVector(vector<double>& list, string delim){
    try {
        string result = "";
        
        if (list.size() == 0) { return result; }
        
        result = list[0];
        
        for (int i = 1; i < list.size(); i++) {
            if (control_pressed) { break;  }
            string temp = toString(list[i]);
            result += delim + temp;
        }
        
        return result;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "getStringFromVector");
        exit(1);
    }
}
/**************************************************************************************************/



