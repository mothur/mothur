/*
 *  mothurOut.cpp
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"
#include "ordervector.hpp"
#include "sharedordervector.h"
#include "counttable.h"

/******************************************************/
MothurOut* MothurOut::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new MothurOut();
	}
	return _uniqueInstance;
}
/*********************************************************************************************/
void MothurOut::appendLogBuffer(string partialLog)  {
    try {
        buffer += partialLog;
    }
    catch(exception& e) {
        errorOut(e, "MothurOut", "appendLogBuffer");
        exit(1);
    }
}
/*********************************************************************************************/
void MothurOut::setLogFileName(string filename, bool append)  {
	try {
        logFileName = filename;
        silenceWarnings = false;
        Utils util;
        if ((filename == "silent")) { silenceLog = true; }
        else {
            if (out.is_open()) { closeLog(); }
            silenceLog = false;
            if (append)     {
                util.openOutputFileAppend(filename, out);
                out << "\n\n************************************************************\n\n\n";
            }else            {  bool opendLog = util.openOutputFile(filename, out);       if (!opendLog) { control_pressed = true; } }
        }
        
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setFileName");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::closeLog()  {
	try {
        if (buffer != "") { string output = buffer; buffer = ""; mothurOut(output);   }
        if (numErrors != 0) {
            if (!silenceLog) {
                out << "\n\n************************************************************\n";
                out << "************************************************************\n";
                out << "************************************************************\n";
                out << "Detected " + toString(numErrors) + " [ERROR] messages, please review.\n";
                out << "************************************************************\n";
                out << "************************************************************\n";
                out << "************************************************************\n";
            }
            logger() << "\n\n************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "Detected " + toString(numErrors) + " [ERROR] messages, please review.\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
            logger() << "************************************************************\n";
        }
        
        if (numWarnings != 0) {
            if (!silenceLog) {
                out << "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                out << "Detected " + toString(numWarnings) + " [WARNING] messages, please review.\n";
                out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
                out << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<^>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
            }
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
        if (buffer != "") { output = buffer + output; buffer = ""; }
        if (output.find("[ERROR]") != string::npos) {
            numErrors++;
            numCommandErrors++;
            if (numCommandErrors > maxCommandErrors) { logger() << "\n**** Exceeded maximum allowed command errors, quitting ****\n"; control_pressed = true; } //abort command
        }
        bool savedSilenceLog = silenceLog;
        bool containsWarning = false;
        if (output.find("[WARNING]") != string::npos) {
            numWarnings++;
            numCommandWarnings++;
            containsWarning = true;
            if (numCommandWarnings > maxCommandWarnings) {
                if (!silenceWarnings) {
                    logger() << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                    out << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                }
                silenceWarnings = true; // write to cout, don't add to logfile
            }
        }
        
        if (!quietMode) {
            if (!silenceLog) {
                if (silenceWarnings && containsWarning) {} //do not print warning to logfile if warnings are silenced
                else { out << output;  }
            }
            logger() << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                if (!silenceLog) { out << output; }
                logger() << output;
            }
        }
        silenceLog = savedSilenceLog;
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToScreen(string output) {
	try {
        if (buffer != "") { output = buffer + output; buffer = ""; }
		if (output.find("[ERROR]") != string::npos) {
            numErrors++;
            numCommandErrors++;
            if (numCommandErrors > maxCommandErrors) { logger() << "\n**** Exceeded maximum allowed command errors, quitting ****\n"; control_pressed = true; } //abort command
        }
        
        bool containsWarning = false;
        if (output.find("[WARNING]") != string::npos) {
            numWarnings++;
            numCommandWarnings++;
            containsWarning = true;
            if (numCommandWarnings > maxCommandWarnings) {
                if (!silenceWarnings) {  logger() << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n"; }
                silenceWarnings = true; // write to cout, don't add to logfile
            }
        }
        
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
            if (!silenceLog) { out << buffer << endl; }
            logger() << buffer << endl;
        }
        buffer = "";
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOut(string output, ofstream& outputFile) {
	try {
        if (buffer != "") { output = buffer + output; buffer = ""; }
        if (output.find("[ERROR]") != string::npos) {
            numErrors++;
            numCommandErrors++;
            if (numCommandErrors > maxCommandErrors) { logger() << "\n**** Exceeded maximum allowed command errors, quitting ****\n"; control_pressed = true; } //abort command
        }
        bool savedSilenceLog = silenceLog;
        bool containsWarning = false;
        if (output.find("[WARNING]") != string::npos) {
            numWarnings++;
            numCommandWarnings++;
            containsWarning = true;
            if (numCommandWarnings > maxCommandWarnings) {
                if (!silenceWarnings) {
                    logger() << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                    out << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                }
                silenceWarnings = true; // write to cout, don't add to logfile
            }
        }
        
        if (!quietMode) {
            if (!silenceLog) {
                if (silenceWarnings && containsWarning) {} //do not print warning to logfile if warnings are silenced
                else { out << output; outputFile << output;  }
            }
            logger() << output;
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                if (!silenceLog) { out << output; }
                outputFile << output;
                logger() << output;
            }
            
        }
        silenceLog = savedSilenceLog;
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
            if (!silenceLog) { out << buffer << endl; }
            logger() << buffer << endl;
            outputFile << buffer << endl;
        }
        buffer = "";
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToLog(string output) {
	try {
        if (buffer != "") { output = buffer + output; buffer = ""; }
        if (output.find("[ERROR]") != string::npos) {
            numErrors++;
            numCommandErrors++;
            if (numCommandErrors > maxCommandErrors) { logger() << "\n**** Exceeded maximum allowed command errors, quitting ****\n"; control_pressed = true; } //abort command
        }
        
        bool savedSilenceLog = silenceLog;
        bool containsWarning = false;
        if (output.find("[WARNING]") != string::npos) {
            numWarnings++;
            numCommandWarnings++;
            containsWarning = true;
            if (numCommandWarnings > maxCommandWarnings) {
                if (!silenceWarnings) {
                    logger() << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                    out << "\n**** Exceeded maximum allowed command warnings, silencing warnings ****\n";
                }
                silenceWarnings = true; // write to cout, don't add to logfile
            }
        }
        
        if (!quietMode) {
            if (!silenceLog) {
                if (silenceWarnings && containsWarning) {} //do not print warning to logfile if warnings are silenced
                else { out << output; }
            }
        }else {
            //check for this being an error
            if ((output.find("[ERROR]") != string::npos) || (output.find("mothur >") != string::npos)) {
                if (!silenceLog) { out << output; }
            }
        }
        silenceLog = savedSilenceLog;
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
    mothurOut("[ERROR]: " + errorType);
    
    unsigned long long ramUsed, total;
    Utils util;
    ramUsed = util.getRAMUsed(); total = util.getTotalRAM();
    mothurOut("RAM used: " + toString(ramUsed/(double)GIG) + "Gigabytes . Total Ram: " + toString(total/(double)GIG) + "Gigabytes.\n\n");
    
    if (pos == string::npos) { //not bad_alloc
        mothurOut(" has occurred in the " + object + " class function " + function + ". Please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry\n");
    }else { //bad alloc
        if (object == "cluster"){
            mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory.  There are two common causes for this, file size and format.\n\nFile Size:\nThe cluster command loads your distance matrix into RAM, and your distance file is most likely too large to fit in RAM. There are two options to help with this. The first is to use a cutoff. By using a cutoff mothur will only load distances that are below the cutoff. If that is still not enough, there is a command called cluster.split, http://www.mothur.org/wiki/cluster.split which divides the distance matrix, and clusters the smaller pieces separately. You may also be able to reduce the size of the original distance matrix by using the commands outlined in the Schloss SOP, http://www.mothur.org/wiki/Schloss_SOP. \n\nWrong Format:\nThis error can be caused by trying to read a column formatted distance matrix using the phylip parameter. By default, the dist.seqs command generates a column formatted distance matrix. To make a phylip formatted matrix set the dist.seqs command parameter output to lt.  \n\nIf you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
        }else if (object == "shhh.flows"){
                mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory. The shhh.flows command is very memory intensive. This error is most commonly caused by trying to process a dataset too large, using multiple processors, or failing to run trim.flows before shhh.flows. If you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required. Running trim.flows with an oligos file, and then shhh.flows with the file option may also resolve the issue. If for some reason you are unable to run shhh.flows with your data, a good alternative is to use the trim.seqs command using a 50-bp sliding window and to trim the sequence when the average quality score over that window drops below 35. Our results suggest that the sequencing error rates by this method are very good, but not quite as good as by shhh.flows and that the resulting sequences tend to be a bit shorter. If you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry. ");
        }else {
            mothurOut(" has occurred in the " + object + " class function " + function + ". This error indicates your computer is running out of memory.  This is most commonly caused by trying to process a dataset too large, using multiple processors, or a file format issue. If you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required. Also, you may be able to reduce the size of your dataset by using the commands outlined in the Schloss SOP, http://www.mothur.org/wiki/Schloss_SOP. If you are unable to resolve the issue, please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
        }
    }
}
/********************************************************************/


