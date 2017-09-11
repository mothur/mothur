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
		void setFileName(string);
		
		void mothurOut(string); //writes to cout and the logfile
		void mothurOutEndLine(); //writes to cout and the logfile
		void mothurOut(string, ofstream&); //writes to the ofstream, cout and the logfile
		void mothurOutEndLine(ofstream&); //writes to the ofstream, cout and the logfile
        void mothurOutJustToScreen(string); //writes to cout
		void mothurOutJustToLog(string);
		void errorOut(exception&, string, string);
		void closeLog();
		string getDefaultPath() { return defaultPath; }
		void setDefaultPath(string);
        string getTestFilePath() { return testFilePath; }
        void setTestFilePath(string);
        string getBlastPath() { return blastPath; }
        void setBlastPath(string);
		string getOutputDir() { return outputDir; }
		void setOutputDir(string);
		
		string getReleaseDate() { return releaseDate; }
		void setReleaseDate(string r) { releaseDate = r; }
		string getVersion() { return version; }
		void setVersion(string r) { version = r; }
	
		void addGroup(string g) { Groups.push_back(g); }
		void setGroups(vector<string>& g) { sort(g.begin(), g.end()); Groups = g; }
		void clearGroups() { Groups.clear(); }
	    int getNumGroups() { return Groups.size(); }
		vector<string> getGroups() { sort(Groups.begin(), Groups.end()); return Groups; }
		void addAllGroup(string g) { namesOfGroups.push_back(g); }
		void setAllGroups(vector<string>& g) { sort(g.begin(), g.end()); namesOfGroups = g; }
		void clearAllGroups() { namesOfGroups.clear(); }
		int getNumAllGroups() { return namesOfGroups.size(); }
	
		vector<string> getAllGroups() { sort(namesOfGroups.begin(), namesOfGroups.end()); return namesOfGroups; }
		
		//functions from mothur.h
		//file operations
        bool dirCheck(string&); //completes path, appends appropriate / or \, makes sure dir is writable.
        bool dirCheck(string&, string); //completes path, appends appropriate / or \, makes sure dir is writable. - no error
        bool mkDir(string&); //completes path, appends appropriate / or \. //returns true it exits or if we can make it
		vector<unsigned long long> divideFile(string, int&); //divides splitting unevenness by sequence
        vector<unsigned long long> divideFile(string filename, int& proc, char delimChar);
        vector<unsigned long long> divideFilePerLine(string, int&); //divides splitting unevenness at line breaks
		int divideFile(string, int&, vector<string>&);
		vector<unsigned long long> setFilePosEachLine(string, int&);
        vector<unsigned long long> setFilePosEachLine(string, unsigned long long&);
		vector<unsigned long long> setFilePosFasta(string, long long&);
        vector<unsigned long long> setFilePosFasta(string, long long&, char);
        vector<unsigned long long> setFilePosFasta(string, int&);
		string sortFile(string, string);
		int appendFiles(string, string);
        int appendBinaryFiles(string, string);
        int appendSFFFiles(string, string);
        int appendFilesWithoutHeaders(string, string);
		int renameFile(string, string); //oldname, newname
		string getFullPathName(string);
        string findProgramPath(string programName);
		string hasPath(string);
		string getExtension(string);
		string getPathName(string);
		string getSimpleName(string);
		string getRootName(string);
		bool isBlank(string);
        bool isCountFile(string);
		int openOutputFile(string, ofstream&);
        int openOutputFileBinary(string, ofstream&);
		int openOutputFileAppend(string, ofstream&);
        int openOutputFileBinaryAppend(string, ofstream&);
		int openInputFile(string, ifstream&);
        int openInputFileBinary(string, ifstream&);
        int openInputFileBinary(string, ifstream&, string);
    #ifdef USE_BOOST
        int openInputFileBinary(string, ifstream&, boost::iostreams::filtering_istream&);
        int openInputFileBinary(string, ifstream&, boost::iostreams::filtering_istream&, string);
    #endif
		int openInputFile(string, ifstream&, string); //no error given
        vector<bool> allGZFiles(vector<string>&);
        vector<bool> isGZ(string); //checks existence and format - will fail for either or both.
        bool fileExists(string name);
        bool checkLocations(string&, string);  //filename, inputDir. checks for file in ./, inputdir, default and mothur's exe location.  Returns false if cant be found. If found completes name with location
		string getline(ifstream&);
		string getline(istringstream&);
        bool stringBlank (string);
		void gobble(istream&);
		void gobble(istringstream&);
        void zapGremlins(istream&);
        void zapGremlins(istringstream&);
        vector<string> splitWhiteSpace(string& rest, char[], int);
        vector<string> splitWhiteSpace(string);
        int splitWhiteSpace(string, vector<float>&, int);
        set<string> readAccnos(string);
        int readAccnos(string, vector<string>&);
        int readAccnos(string, vector<string>&, string);
        map<string, int> readNames(string);
        map<string, int> readNames(string, unsigned long int&);
        int readTax(string, map<string, string>&, bool);
        vector<consTax> readConsTax(string);
        int readConsTax(string, map<int, consTax2>&);
        int readNames(string, map<string, string>&, map<string, int>&);
		int readNames(string, map<string, string>&);
        int readNames(string, map<string, string>&, bool);
        int readNames(string, map<string, string>&, int);
		int readNames(string, map<string, vector<string> >&);
		int readNames(string, vector<seqPriorityNode>&, map<string, string>&);
		int mothurRemove(string);
        int printVsearchFile(vector<seqPriorityNode>&, string, string, string); //sorts and prints by abundance adding /ab=xxx/
        bool mothurConvert(char, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
		bool mothurConvert(string, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
        bool mothurConvert(string, intDist&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
		bool mothurConvert(string, float&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
		bool mothurConvert(string, double&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
        bool mothurConvert(char, string&);
		
		//searchs and checks
		bool checkReleaseVersion(ifstream&, string);
        int getTimeStamp(string filename);
		bool anyLabelsToProcess(string, set<string>&, string);
		bool inUsersGroups(vector<string>, vector<string>); //returns true if any of the strings in first vector are in second vector
        bool inUsersGroups(vector<int>, vector< vector<int> >);
		bool inUsersGroups(string, vector<string>);
        bool inUsersGroups(int, vector<int>);
		void getNumSeqs(ifstream&, int&);
		int getNumSeqs(ifstream&);
		int getNumNames(string);
		int getNumChar(string, char);
		bool isTrue(string);
		bool isContainingOnlyDigits(string);
        bool containsAlphas(string);
        bool isAllAlphas(string);
        bool isAllAlphaNumerics(string);
		bool isNumeric1(string);
        bool isNumeric1(char);
        bool isInteger(string);
        bool allSpaces(string);
        bool isLabelEquivalent(string, string);
        string getSimpleLabel(string);
        string findEdianness();
        string mothurGetpid(int);
        unsigned long long getRAMUsed();
        unsigned long long getTotalRAM();
        string getStringFromVector(vector<string>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
        string getStringFromVector(vector<int>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
        string getStringFromVector(vector<double>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
		
		//string manipulation
		void splitAtEquals(string&, string&);
		void splitAtComma(string&, string&);	
		void splitAtComma(string&, vector<string>&);
		void splitAtDash(string&, set<int>&);
		void splitAtDash(string&, set<string>&);
		void splitAtDash(string&, vector<string>&);
		void splitAtChar(string&, vector<string>&, char);
        void splitAtChar(string&, string&, char);
        int removeBlanks(vector<string>&);
        vector<string> splitWhiteSpaceWithQuotes(string);
        int splitWhiteSpaceWithQuotes(string, vector<string>&);
		float removeConfidences(string&);
        string removeQuotes(string);
        string makeList(vector<string>&);
        bool isSubset(vector<string>, vector<string>); //bigSet, subset
        int checkName(string&);
        bool checkGroupName(string name);
        map<string, vector<string> > parseClasses(string);
        string addUnclassifieds(string tax, int maxlevel, bool probs);
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
    
		
		//math operation
        int max(int, int);
        int min(int, int);
        double max(vector<double>&); //returns largest value in vector
        double min(vector<double>&); //returns smallest value in vector
        int max(vector<int>);
        float max(vector<float>);
        int min(vector<int>);
        float min(vector<float>);
		int factorial(int num);
		vector<vector<double> > binomial(int);
		float ceilDist(float, int);
		float roundDist(float, int);
		unsigned int fromBase36(string);
        double median(vector<double>);
        int median(vector<int>);
        int average(vector<int>);
        int sum(vector<int>);
        double sum(vector<double>);
    
		int getRandomIndex(int); //highest
        int getRandomNumber();
        double getRandomDouble0to1();
        double getStandardDeviation(vector<int>&);
        vector<double> getStandardDeviation(vector< vector<double> >&);
        vector<double> getStandardDeviation(vector< vector<double> >&, vector<double>&);
        vector<double> getAverages(vector< vector<double> >&);
        double getAverage(vector<double>);
        vector< vector<seqDist> > getStandardDeviation(vector< vector< vector<seqDist> > >&);
        vector< vector<seqDist> > getStandardDeviation(vector< vector< vector<seqDist> > >&, vector< vector<seqDist> >&);
        vector< vector<seqDist> > getAverages(vector< vector< vector<seqDist> > >&, string);
        vector< vector<seqDist> > getAverages(vector< vector< vector<seqDist> > >&);
    
        int getControl_pressed()                        { return control_pressed;           }
        void setControl_pressed(int t)                  { control_pressed = t;              }
        vector<string> getTreenames()                   { return Treenames;                 }
        void setTreenames(vector<string> t)             { Treenames = t;                    }
        vector<string> getSharedBinLabelsInFile()       { return sharedBinLabelsInFile;     }
        void setSharedBinLabelsInFile(vector<string> t) { sharedBinLabelsInFile = t;        }
        vector<string> getCurrentSharedBinLabels()      { return currentSharedBinLabels;    }
        void setCurrentSharedBinLabels(vector<string> t){ currentSharedBinLabels = t;       }
        void eraseCurrentSharedBinLabel(int bin)        { currentSharedBinLabels.erase(currentSharedBinLabels.begin()+bin);       }
        vector<string> getListBinLabelsInFile()         { return listBinLabelsInFile;       }
        void setListBinLabelsInFile(vector<string> t)   { listBinLabelsInFile = t;          }
        string getSaveNextLabel()                       { return saveNextLabel;             }
        void setSaveNextLabel(string t)                 { saveNextLabel = t;                }
        string getSharedHeaderMode()                    { return sharedHeaderMode;          }
        void setSharedHeaderMode(string t)              { sharedHeaderMode = t;             }
        string getGroupMode()                           { return groupMode;                 }
        void setGroupMode(string t)                     { groupMode = t;                    }
        string getTestDirectory()                       { return testDirectory;             }
        void setTestDirectory(string t)                 { testDirectory = t;                }
        bool getPrintedSharedHeaders()                  { return printedSharedHeaders;      }
        void setPrintedSharedHeaders(bool t)            { printedSharedHeaders = t;         }
        bool getPrintedListHeaders()                    { return printedListHeaders;        }
        void setPrintedListHeaders(bool t)              { printedListHeaders = t;           }
        bool getCommandInputsConvertError()             { return commandInputsConvertError; }
        void setCommandInputsConvertError(bool t)       { commandInputsConvertError = t;    }
        bool getChangedSeqNames()                       { return changedSeqNames;           }
        void setChangedSeqNames(bool t)                 { changedSeqNames = t;              }
        bool getModifyNames()                           { return modifyNames;               }
        void setModifyNames(bool t)                     { modifyNames = t;                  }
        bool getExecuting()                             { return executing;                 }
        void setExecuting(bool t)                       { executing = t;                    }
        bool getRunParse()                              { return runParse;                  }
        void setRunParse(bool t)                        { runParse = t;                     }
        bool getJumble()                                { return jumble;                    }
        void setJumble(bool t)                          { jumble = t;                       }
        bool getGui()                                   { return gui;                       }
        void setGui(bool t)                             { gui = t;                          }
        bool getMothurCalling()                         { return mothurCalling;             }
        void setMothurCalling(bool t)                   { mothurCalling = t;                }
        bool getDebug()                                 { return debug;                     }
        void setDebug(bool t)                           { debug = t;                        }
        bool getQuietMode()                             { return quietMode;                 }
        void setQuietMode(bool t)                       { quietMode = t;                    }

    
		//current files - if you add a new type you must edit optionParser->getParameters, get.current and set.current commands and mothurOut->printCurrentFiles/clearCurrentFiles/getCurrentTypes/hasCurrentFiles. add a get and set function.
		string getPhylipFile()		{ return phylipfile;		}
		string getColumnFile()		{ return columnfile;		}
		string getListFile()		{ return listfile;			}
		string getRabundFile()		{ return rabundfile;		}
		string getSabundFile()		{ return sabundfile;		}
		string getNameFile()		{ return namefile;			}	
		string getGroupFile()		{ return groupfile;			}	
		string getOrderFile()		{ return orderfile;			}
		string getOrderGroupFile()	{ return ordergroupfile;	}
		string getTreeFile()		{ return treefile;			}
		string getSharedFile()		{ return sharedfile;		}
		string getRelAbundFile()	{ return relabundfile;		}
		string getDesignFile()		{ return designfile;		}
		string getFastaFile()		{ return fastafile;			}
		string getSFFFile()			{ return sfffile;			}
		string getQualFile()		{ return qualfile;			}
		string getOligosFile()		{ return oligosfile;		}
		string getAccnosFile()		{ return accnosfile;		}
		string getTaxonomyFile()	{ return taxonomyfile;		}
		string getFlowFile()		{ return flowfile;			}
        string getBiomFile()		{ return biomfile;			}
        string getCountTableFile()	{ return counttablefile;	}
        string getSummaryFile()     { return summaryfile;       }
        string getFileFile()        { return filefile;          }
		string getProcessors()		{ return processors;		}
        string getConsTaxonomyFile(){ return constaxonomyfile;  }
        string getContigsReportFile(){ return contigsreportfile;  }
        int getNumErrors()          { return numErrors;         }
        string getProgramPath()     { return mothurProgramPath;  }
		
		void setListFile(string f)			{ listfile = getFullPathName(f);			}
		void setTreeFile(string f)			{ treefile = getFullPathName(f);			}
		void setGroupFile(string f)			{ groupfile = getFullPathName(f);	groupMode = "group";		}		
		void setPhylipFile(string f)		{ phylipfile = getFullPathName(f);			}
		void setColumnFile(string f)		{ columnfile = getFullPathName(f);			}
		void setNameFile(string f)			{ namefile = getFullPathName(f);			}	
		void setRabundFile(string f)		{ rabundfile = getFullPathName(f);			}
		void setSabundFile(string f)		{ sabundfile = getFullPathName(f);			}
		void setSharedFile(string f)		{ sharedfile = getFullPathName(f);			}
		void setRelAbundFile(string f)		{ relabundfile = getFullPathName(f);		}
		void setOrderFile(string f)			{ orderfile = getFullPathName(f);			}
		void setOrderGroupFile(string f)	{ ordergroupfile = getFullPathName(f);		}
		void setDesignFile(string f)		{ designfile = getFullPathName(f);			}
		void setFastaFile(string f)			{ fastafile = getFullPathName(f);			}
		void setSFFFile(string f)			{ sfffile = getFullPathName(f);				}
		void setQualFile(string f)			{ qualfile = getFullPathName(f);			}
		void setOligosFile(string f)		{ oligosfile = getFullPathName(f);			}
		void setAccnosFile(string f)		{ accnosfile = getFullPathName(f);			}
		void setTaxonomyFile(string f)		{ taxonomyfile = getFullPathName(f);		}
		void setFlowFile(string f)			{ flowfile = getFullPathName(f);			}
        void setBiomFile(string f)			{ biomfile = getFullPathName(f);			}
        void setSummaryFile(string f)		{ summaryfile = getFullPathName(f);			}
        void setFileFile(string f)          { filefile = getFullPathName(f);			}
        void setConsTaxonomyFile(string f)  { constaxonomyfile = getFullPathName(f);	}
        void setContigsReportFile(string f) { contigsreportfile = getFullPathName(f);	}
        void setCountTableFile(string f)	{ counttablefile = getFullPathName(f);	groupMode = "count";	}
        void setProcessors(string p)		{ processors = p; mothurOut("\nUsing " + toString(p) + " processors.\n");	}
        void setProgramPath(string);
		
		void printCurrentFiles(string); //string="" for just to logfile.
		bool hasCurrentFiles();
		void clearCurrentFiles();
        set<string> getCurrentTypes(); 
		
	private:
		static MothurOut* _uniqueInstance;
		MothurOut( const MothurOut& ); // Disable copy constructor
		void operator=( const MothurOut& ); // Disable assignment operator
		MothurOut() { 
			control_pressed = false; defaultPath=""; blastPath=""; testFilePath = "";
            filefile = "";
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
            constaxonomyfile = "";
			processors = "1";
			flowfile = "";
            biomfile = "";
            counttablefile = "";
            summaryfile = "";
            contigsreportfile = "";
			gui = false;
			printedSharedHeaders = false;
            printedListHeaders = false;
			commandInputsConvertError = false;
            mothurCalling = false;
            debug = false;
            quietMode = false;
			sharedHeaderMode = "";
            groupMode = "group";
            changedSeqNames = false;
            modifyNames = true;
            numErrors = 0;
            numWarnings = 0;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            mersenne_twister_engine.seed(seed);
		}
		~MothurOut();

		string logFileName, mothurProgramPath;
		string defaultPath, outputDir, blastPath;
		string releaseDate, version;
	
		string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile, biomfile, filefile, testFilePath;
		string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, processors, flowfile, counttablefile, summaryfile, constaxonomyfile, contigsreportfile;
        mt19937_64 mersenne_twister_engine;

		vector<string> Groups;
		vector<string> namesOfGroups;
		ofstream out;
        int numErrors, numWarnings;
        mutex threadModerator;
        vector<string> Treenames;
        vector<string> sharedBinLabelsInFile;
        vector<string> currentSharedBinLabels;
        vector<string> listBinLabelsInFile;
        string saveNextLabel, sharedHeaderMode, groupMode, testDirectory;
        bool printedSharedHeaders, printedListHeaders, commandInputsConvertError, changedSeqNames, modifyNames;
        bool executing, runParse, jumble, gui, mothurCalling, debug, quietMode;
        int control_pressed;
		
};
/***********************************************/

#endif

