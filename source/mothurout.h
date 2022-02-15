#ifndef MOTHUROUT_H
#define MOTHUROUT_H

/*
 *  mothurOut.h
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"

/***********************************************/
struct logger {
    
    logger() {}
    ~logger() {}
    
    template< class T >
    logger& operator <<( const T& o ) {
        //lock_guard<std::mutex> guard(token);
        cout << o; return *this;
    }
    
    logger& operator<<(ostream& (*m)(ostream&) ) {
        //lock_guard<std::mutex> guard(token);
        cout << m; return *this;
    }
private:
    //std::mutex token;
};
/***********************************************/

class MothurOut {
	
	public:
		static MothurOut* getInstance();
    
        //logger
        bool getDebug()                                 { return debug;                                 }
        void setDebug(bool t)                           { debug = t;                                    }
        void setQuietMode(bool t)                       { quietMode = t;                                }
        int getNumErrors()                              { return numErrors;                             }
        void resetCommandErrors()                       { control_pressed = false; numCommandErrors = 0; numCommandWarnings = 0; silenceWarnings = false;}
        string getLogFileName()                         { return logFileName;                           }
		void setLogFileName(string f, bool append);
        void setHomePath(string);
        string getHomePath()  { return homePath; }
        void setPaths(vector<string>); //environment variable 'PATH' values
        vector<string> getPaths() { return paths; }
    
		void mothurOut(string); //writes to cout and the logfile
		void mothurOutEndLine(); //writes to cout and the logfile
		void mothurOut(string, ofstream&); //writes to the ofstream, cout and the logfile
        void mothurOutJustToScreen(string); //writes to cout
		void mothurOutJustToLog(string);
		void errorOut(exception&, string, string);
		void closeLog();
    
        //globals
        void setRandomSeed(unsigned s)                  { seed = s;                         }
        unsigned getRandomSeed()                        { return seed;                      }
        bool getControl_pressed()                       { return control_pressed;           }
        void setControl_pressed(bool t)                 { control_pressed = t;              }
        bool getChangedSeqNames()                       { return changedSeqNames;           }
        void setChangedSeqNames(bool t)                 { changedSeqNames = t;              }
        bool getChangedGroupNames()                     { return changedGroupNames;         }
        void setChangedGroupNames(bool t)               { changedGroupNames = t;            }
        bool getExecuting()                             { return executing;                 }
        void setExecuting(bool t)                       { executing = t;                    }
    
        vector<vector<vector<char>>> codons;
        set<char> validAminoAcids;
    
	private:
		static MothurOut* _uniqueInstance;
		MothurOut( const MothurOut& ); // Disable copy constructor
		void operator=( const MothurOut& ); // Disable assignment operator
		MothurOut() { 
			control_pressed = false;
            debug = false;
            quietMode = false;
            changedSeqNames = true;
            changedGroupNames = true;
            silenceLog = false;
            silenceWarnings = false; //per command
            numErrors = 0; numWarnings = 0;
            numCommandErrors = 0; numCommandWarnings = 0;
            maxCommandErrors = 10; maxCommandWarnings = 10;
            logFileName = "";
            buffer = "";
            homePath = "";
            outLog = NULL;
            seed = std::chrono::system_clock::now().time_since_epoch().count();
            
            initialize(); //fills validAminoAcids and codons
		}
		~MothurOut();
		
        void appendLogBuffer(string); //used to store log before we establish the logfilename
        void initialize();

		ofstream* outLog;
        unsigned seed;
        int numErrors, numWarnings, numCommandErrors, numCommandWarnings, maxCommandErrors, maxCommandWarnings;
        string logFileName, buffer, homePath;
        vector<string> paths;
        bool changedSeqNames, changedGroupNames, silenceLog, silenceWarnings, control_pressed, executing, debug, quietMode;
};
/***********************************************/


#endif

