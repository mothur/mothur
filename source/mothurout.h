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
		vector<string> Treenames;
		vector<string> sharedBinLabelsInFile;
		vector<string> currentSharedBinLabels;
        vector<string> listBinLabelsInFile;
		string saveNextLabel, argv, sharedHeaderMode, groupMode;
		bool printedSharedHeaders, printedListHeaders, commandInputsConvertError, changedSeqNames, modifyNames;
		
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
		bool mothurConvert(string, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
        bool mothurConvert(string, intDist&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
		bool mothurConvert(string, float&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
		bool mothurConvert(string, double&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
	
		
		//searchs and checks
		bool checkReleaseVersion(ifstream&, string);
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
		bool isNumeric1(string);
        bool isInteger(string);
        bool isLabelEquivalent(string, string);
        string getSimpleLabel(string);
        string findEdianness();
        string mothurGetpid(int);
	
		
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
		int removeConfidences(string&);
        string removeQuotes(string);
        string makeList(vector<string>&);
        bool isSubset(vector<string>, vector<string>); //bigSet, subset
        int checkName(string&);
        map<string, vector<string> > parseClasses(string);
        
		
		//math operation
        int max(int, int);
        int min(int, int);
        double max(vector<double>&); //returns largest value in vector
        double min(vector<double>&); //returns smallest value in vector
		int factorial(int num);
		vector<vector<double> > binomial(int);
		float ceilDist(float, int);
		float roundDist(float, int);
		unsigned int fromBase36(string);
        double median(vector<double>);
        int median(vector<int>);
        int average(vector<int>);
        int sum(vector<int>);
		int getRandomIndex(int); //highest
        double getStandardDeviation(vector<int>&);
        vector<double> getStandardDeviation(vector< vector<double> >&);
        vector<double> getStandardDeviation(vector< vector<double> >&, vector<double>&);
        vector<double> getAverages(vector< vector<double> >&);
        double getAverage(vector<double>);
        vector< vector<seqDist> > getStandardDeviation(vector< vector< vector<seqDist> > >&);
        vector< vector<seqDist> > getStandardDeviation(vector< vector< vector<seqDist> > >&, vector< vector<seqDist> >&);
        vector< vector<seqDist> > getAverages(vector< vector< vector<seqDist> > >&, string);
        vector< vector<seqDist> > getAverages(vector< vector< vector<seqDist> > >&);

		int control_pressed;
		bool executing, runParse, jumble, gui, mothurCalling, debug, quietMode;
		
		//current files - if you add a new type you must edit optionParser->getParameters, get.current and set.current commands and mothurOut->printCurrentFiles/clearCurrentFiles/getCurrentTypes. add a get and set function.
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
        void setCountTableFile(string f)	{ counttablefile = getFullPathName(f);	groupMode = "count";	}
        void setProcessors(string p)		{ processors = p; mothurOut("\nUsing " + toString(p) + " processors.\n");	}
		
		void printCurrentFiles();
		bool hasCurrentFiles();
		void clearCurrentFiles();
        set<string> getCurrentTypes(); 
		
	private:
		static MothurOut* _uniqueInstance;
		MothurOut( const MothurOut& ); // Disable copy constructor
		void operator=( const MothurOut& ); // Disable assignment operator
		MothurOut() { 
			control_pressed = false; defaultPath="";
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
			processors = "1";
			flowfile = "";
            biomfile = "";
            counttablefile = "";
            summaryfile = "";
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
		}
		~MothurOut();

		string logFileName;
		string defaultPath, outputDir;
		string releaseDate, version;
	
		string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile, biomfile, filefile;
		string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, processors, flowfile, counttablefile, summaryfile;

		vector<string> Groups;
		vector<string> namesOfGroups;
		ofstream out;
		
		int mem_usage(double&, double&);

};
/***********************************************/

#endif

