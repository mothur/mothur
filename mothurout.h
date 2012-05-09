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
		//map<string, string> names;
		vector<string> binLabelsInFile;
		vector<string> currentBinLabels;
		string saveNextLabel, argv, sharedHeaderMode;
		bool printedHeaders, commandInputsConvertError;
		
		//functions from mothur.h
		//file operations
        bool dirCheck(string&); //completes path, appends appropriate / or \, makes sure dir is writable.
		vector<unsigned long long> divideFile(string, int&);
		int divideFile(string, int&, vector<string>&);
		vector<unsigned long long> setFilePosEachLine(string, int&);
		vector<unsigned long long> setFilePosFasta(string, int&);
		string sortFile(string, string);
		int appendFiles(string, string);
		int renameFile(string, string); //oldname, newname
		string getFullPathName(string);
		string hasPath(string);
		string getExtension(string);
		string getPathName(string);
		string getSimpleName(string);
		string getRootName(string);
		bool isBlank(string);
		int openOutputFile(string, ofstream&);
		int openOutputFileAppend(string, ofstream&);
		int openInputFile(string, ifstream&);
		int openInputFile(string, ifstream&, string); //no error given 
		string getline(ifstream&);
		string getline(istringstream&);
		void gobble(istream&);
		void gobble(istringstream&);
		map<string, int> readNames(string);
		int readNames(string, map<string, string>&);
		int readNames(string, map<string, vector<string> >&);
		int readNames(string, vector<seqPriorityNode>&, map<string, string>&);
		int mothurRemove(string);
		bool mothurConvert(string, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
		bool mothurConvert(string, float&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
		bool mothurConvert(string, double&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
	
		
		//searchs and checks
		bool checkReleaseVersion(ifstream&, string);
		bool anyLabelsToProcess(string, set<string>&, string);
		bool inUsersGroups(vector<string>, vector<string>);
		bool inUsersGroups(string, vector<string>);
		void getNumSeqs(ifstream&, int&);
		int getNumSeqs(ifstream&);
		int getNumNames(string);
		int getNumChar(string, char);
		bool isTrue(string);
		bool isContainingOnlyDigits(string);
		bool isNumeric1(string);
	
		
		//string manipulation
		void splitAtEquals(string&, string&);
		void splitAtComma(string&, string&);	
		void splitAtComma(string&, vector<string>&);
		void splitAtDash(string&, set<int>&);
		void splitAtDash(string&, set<string>&);
		void splitAtDash(string&, vector<string>&);
		void splitAtChar(string&, vector<string>&, char);
        void splitAtChar(string&, string&, char);
		int removeConfidences(string&);
		
		//math operation
		int factorial(int num);
		vector<vector<double> > binomial(int);
		float ceilDist(float, int);
		float roundDist(float, int);
		unsigned int fromBase36(string);
		int getRandomIndex(int); //highest

		int control_pressed;
		bool executing, runParse, jumble, gui, mothurCalling, debug;
		
		//current files - if you add a new type you must edit optionParser->getParameters, get.current command and mothurOut->printCurrentFiles/clearCurrentFiles.
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
		string getProcessors()		{ return processors;		}
		
		void setListFile(string f)			{ listfile = getFullPathName(f);			}
		void setTreeFile(string f)			{ treefile = getFullPathName(f);			}
		void setGroupFile(string f)			{ groupfile = getFullPathName(f);			}		
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
		void setProcessors(string p)		{ processors = p;							}
		
		void printCurrentFiles();
		bool hasCurrentFiles();
		void clearCurrentFiles();
		
	private:
		static MothurOut* _uniqueInstance;
		MothurOut( const MothurOut& ); // Disable copy constructor
		void operator=( const MothurOut& ); // Disable assignment operator
		MothurOut() { 
			control_pressed = false; defaultPath=""; 
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
			gui = false;
			printedHeaders = false;
			commandInputsConvertError = false;
            mothurCalling = false;
            debug = false;
			sharedHeaderMode = "";
		}
		~MothurOut();

		string logFileName;
		string defaultPath, outputDir;
		string releaseDate, version;
	
		string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile, biomfile;
		string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, processors, flowfile;

		vector<string> Groups;
		vector<string> namesOfGroups;
		ofstream out;
		
		int mem_usage(double&, double&);

};
/***********************************************/

#endif

