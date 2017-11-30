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
        cout << o; return *this;
    }
    
    logger& operator<<(ostream& (*m)(ostream&) ) {
        cout << m; return *this;
    }
    
}; 
/***********************************************/
class OrderVector;
class SharedOrderVector;
class RAbundVector;
class SharedRAbundVector;

class MothurOut {
	
	public:
		static MothurOut* getInstance();
    
        //logger
        void appendLogBuffer(string); //used to store log before we establish the logfilename
        bool getDebug()                                 { return debug;                     }
        void setDebug(bool t)                           { debug = t;                        }
        bool getQuietMode()                             { return quietMode;                 }
        void setQuietMode(bool t)                       { quietMode = t;                    }
        int getNumErrors()          { return numErrors;         }
		void setLogFileName(string f, bool append);
        string getLogFileName() { return logFileName; }
		void mothurOut(string); //writes to cout and the logfile
		void mothurOutEndLine(); //writes to cout and the logfile
		void mothurOut(string, ofstream&); //writes to the ofstream, cout and the logfile
		void mothurOutEndLine(ofstream&); //writes to the ofstream, cout and the logfile
        void mothurOutJustToScreen(string); //writes to cout
		void mothurOutJustToLog(string);
		void errorOut(exception&, string, string);
		void closeLog();
    
        //random operations
        int getRandomIndex(int); //highest
        int getRandomNumber();
        double getRandomDouble0to1();
        int mothurRandomShuffle(vector<int>&);
        int mothurRandomShuffle(vector< vector<double> >&);
        int mothurRandomShuffle(vector<string>&);
        int mothurRandomShuffle(vector<item>&);
        int mothurRandomShuffle(vector<PCell*>&);
        int mothurRandomShuffle(vector<PDistCellMin>&);
        int mothurRandomShuffle(OrderVector&);
        int mothurRandomShuffle(SharedOrderVector&);
        int mothurRandomShuffle(vector<SharedRAbundVector*>&);
        void setRandomSeed(unsigned s) { mersenne_twister_engine.seed(s); srand(s); }
    
        //globals
        
        bool getControl_pressed()                       { return control_pressed;           }
        void setControl_pressed(bool t)                 { control_pressed = t;              }
        bool getPrintedSharedHeaders()                  { return printedSharedHeaders;      }
        void setPrintedSharedHeaders(bool t)            { printedSharedHeaders = t;         }
        bool getPrintedListHeaders()                    { return printedListHeaders;        }
        void setPrintedListHeaders(bool t)              { printedListHeaders = t;           }
        bool getChangedSeqNames()                       { return changedSeqNames;           }
        void setChangedSeqNames(bool t)                 { changedSeqNames = t;              }
        bool getModifyNames()                           { return modifyNames;               }
        void setModifyNames(bool t)                     { modifyNames = t;                  }
        bool getExecuting()                             { return executing;                 }
        void setExecuting(bool t)                       { executing = t;                    }
    
    bool getJumble()                                { return jumble;                    }
        void setJumble(bool t)                          { jumble = t;                       }
        bool getMothurCalling()                         { return mothurCalling;             }
        void setMothurCalling(bool t)                   { mothurCalling = t;                }
    
    //should be owned by datastructures and passed
   // vector<string> getTreenames()                   { return Treenames;                 }
  //  void setTreenames(vector<string> t)             { Treenames = t;                    }
  //  string getSaveNextLabel()                       { return saveNextLabel;             }
  //  void setSaveNextLabel(string t)                 { saveNextLabel = t;                }
  //      
    
	
				
	private:
		static MothurOut* _uniqueInstance;
		MothurOut( const MothurOut& ); // Disable copy constructor
		void operator=( const MothurOut& ); // Disable assignment operator
		MothurOut() { 
			control_pressed = false;
			printedSharedHeaders = false;
            printedListHeaders = false;
            mothurCalling = false;
            debug = false;
            quietMode = false;
			sharedHeaderMode = "";
            changedSeqNames = false;
            modifyNames = true;
            numErrors = 0;
            numWarnings = 0;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            mersenne_twister_engine.seed(seed);
            logFileName = ""; buffer = "";
		}
		~MothurOut();

		
        mt19937_64 mersenne_twister_engine;
		ofstream out;
        int numErrors, numWarnings;
        //vector<string> Treenames;
        string sharedHeaderMode,  logfileName, logFileName;
        bool printedSharedHeaders, printedListHeaders, changedSeqNames, modifyNames;
        bool executing, runParse, jumble, mothurCalling, debug, quietMode;
        bool control_pressed;
        string buffer;
        //std::mutex token;
    
		
};
/***********************************************/

#endif

