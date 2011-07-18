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
		vector<string> Groups;
		vector<string> Treenames;
		map<string, string> names;
		vector<string> namesOfGroups;
		vector<string> binLabelsInFile;
		vector<string> currentBinLabels;
		string saveNextLabel, argv, sharedHeaderMode;
		bool printedHeaders;
		
		//functions from mothur.h
		//file operations
		vector<unsigned long int> divideFile(string, int&);
		int divideFile(string, int&, vector<string>&);
		vector<unsigned long int> setFilePosEachLine(string, int&);
		vector<unsigned long int> setFilePosFasta(string, int&);
		string sortFile(string, string);
		void appendFiles(string, string);
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
		int readNames(string, vector<seqPriorityNode>&, map<string, string>&);
		void mothurRemove(string);
		
		//searchs and checks
		bool checkReleaseVersion(ifstream&, string);
		bool anyLabelsToProcess(string, set<string>&, string);
		bool inUsersGroups(vector<string>, vector<string>);
		bool inUsersGroups(string, vector<string>);
		void getNumSeqs(ifstream&, int&);
		int getNumSeqs(ifstream&);
		int getNumNames(string);
		bool isTrue(string);
		bool isContainingOnlyDigits(string);
	
		
		//string manipulation
		void splitAtEquals(string&, string&);
		void splitAtComma(string&, string&);	
		void splitAtComma(string&, vector<string>&);
		void splitAtDash(string&, set<int>&);
		void splitAtDash(string&, set<string>&);
		void splitAtDash(string&, vector<string>&);
		void splitAtChar(string&, vector<string>&, char);
		
		//math operation
		int factorial(int num);
		vector<vector<double> > binomial(int);
		float ceilDist(float, int);
		float roundDist(float, int);
		unsigned int fromBase36(string);
		int getRandomIndex(int); //highest

		int control_pressed;
		bool executing, runParse, jumble, gui;
		
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
			gui = false;
			printedHeaders = false;
			sharedHeaderMode = "";
		}
		~MothurOut();

		string logFileName;
		string defaultPath, outputDir;
		string releaseDate, version;
	
		string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile;
		string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, processors, flowfile;

	
		ofstream out;
		
		int mem_usage(double&, double&);

};
/***********************************************/

#endif

