#ifndef CURRENTFILE_H
#define CURRENTFILE_H

/*
 *  currentfile.h
 *  Mothur
 *
 *  Created by westcott on 3/15/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


//NOT ThreadSafe - but designed to be read only from threads and read write from main thread.

#include "mothurout.h"
#include "utils.hpp"

/***********************************************/

class CurrentFile {
	
	public:
		static CurrentFile* getInstance() {
			if(instance == 0) {	instance = new CurrentFile();	}
			return instance;
		}
		
    
        unsigned long long getRAMUsed();
        unsigned long long getTotalRAM();
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
		
		void setListFile(string f)			{ listfile = util.getFullPathName(f);			}
        void setBiomFile(string f)			{ biomfile = util.getFullPathName(f);			}
        void setFlowFile(string f)			{ flowfile = util.getFullPathName(f);			}
        void setContigsReportFile(string f)	{ contigsreportfile = util.getFullPathName(f);			}
        void setSummaryFile(string f)		{ summaryfile = util.getFullPathName(f);		}
		void setTreeFile(string f)			{ treefile = util.getFullPathName(f);			}
        void setGroupFile(string f)			{ groupfile = util.getFullPathName(f);	setGroupMode("group");		}
        void setCountFile(string f)			{ countfile = util.getFullPathName(f);	setGroupMode("count");		}
		void setPhylipFile(string f)		{ phylipfile = util.getFullPathName(f);			}
		void setColumnFile(string f)		{ columnfile = util.getFullPathName(f);			}
		void setNameFile(string f)			{ namefile = util.getFullPathName(f);			}
		void setRabundFile(string f)		{ rabundfile = util.getFullPathName(f);			}
		void setSabundFile(string f)		{ sabundfile = util.getFullPathName(f);			}
		void setSharedFile(string f)		{ sharedfile = util.getFullPathName(f);			}
		void setRelAbundFile(string f)		{ relabundfile = util.getFullPathName(f);		}
		void setOrderFile(string f)			{ orderfile = util.getFullPathName(f);			}
		void setOrderGroupFile(string f)	{ ordergroupfile = util.getFullPathName(f);		}
		void setDesignFile(string f)		{ designfile = util.getFullPathName(f);			}
		void setFastaFile(string f)			{ fastafile = util.getFullPathName(f);			}
		void setSFFFile(string f)			{ sfffile = util.getFullPathName(f);			}
		void setQualFile(string f)			{ qualfile = util.getFullPathName(f);			}
		void setOligosFile(string f)		{ oligosfile = util.getFullPathName(f);			}
        void setAccnosFile(string f)		{ accnosfile = util.getFullPathName(f);			}
        void setTaxonomyFile(string f)		{ taxonomyfile = util.getFullPathName(f);       }
        void setConsTaxonomyFile(string f)  { constaxonomyfile = util.getFullPathName(f);	}
        void setProgramPath(string f)       { mothurProgramPath = util.getFullPathName(f);	}
        void setFileFile(string f)          { filefile = util.getFullPathName(f);           }
    
        //current files - if you add a new type you must edit optionParser->getParameters, get.current and set.current commands and mothurOut->printCurrentFiles/clearCurrentFiles/getCurrentTypes/hasCurrentFiles. add a get and set function.
        string getAccnosFile()		{ return accnosfile;		}
        string getTaxonomyFile()	{ return taxonomyfile;		}
        string getFlowFile()		{ return flowfile;			}
        string getContigsReportFile(){ return contigsreportfile;			}
        string getBiomFile()		{ return biomfile;			}
        string getCountFile()       { return countfile;         }
        string getSummaryFile()     { return summaryfile;       }
        string getFileFile()        { return filefile;          }
        string getProcessors()		{ return processors;		}
        int setProcessors(string p);
        string getConsTaxonomyFile(){ return constaxonomyfile;  }
        string getProgramPath()     { return mothurProgramPath;  }
    
        string getDefaultPath() { return defaultPath; }
        void setDefaultPath(string);
        string getTestFilePath() { return testFilePath; }
        void setTestFilePath(string);
        string getBlastPath() { return blastPath; }
        void setBlastPath(string);
        string getOutputDir() { return outputDir; }
        void setOutputDir(string f) { outputDir = util.getFullPathName(f); }
        string getInputDir() { return inputDir; }
        void setInputDir(string f) { inputDir = util.getFullPathName(f); }
        void setFileName(string);
        string getReleaseDate() { return releaseDate; }
        void setReleaseDate(string r) { releaseDate = r; }
        string getVersion() { return version; }
        void setVersion(string r) { version = r; }
        vector<string> getLocations() { vector<string> locations; locations.push_back(inputDir); locations.push_back(outputDir); locations.push_back(defaultPath); locations.push_back(mothurProgramPath); return locations; }
    


        void printCurrentFiles(string); //string="" for just to logfile.
        void clearCurrentFiles();
        set<string> getCurrentTypes();
        bool hasCurrentFiles();
    
        string getSharedHeaderMode()                    { return sharedHeaderMode;          }
        void setSharedHeaderMode(string t)              { sharedHeaderMode = t;             }
        string getGroupMode()                           { return groupMode;                 }
        void setGroupMode(string t)                     { groupMode = t;                    }
        string getTestDirectory()                       { return testDirectory;             }
        void setTestDirectory(string t)                 { testDirectory = t;                }
   
	
	private:
		MothurOut* m;
        Utils util;
        string logFileName, mothurProgramPath;
        string defaultPath, outputDir, blastPath, inputDir;
        string releaseDate, version;
    
        string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile, biomfile, filefile, testFilePath, contigsreportfile;
        string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, processors, flowfile, countfile, summaryfile, constaxonomyfile, groupMode, testDirectory, sharedHeaderMode;
    
		
		static CurrentFile* instance;
		CurrentFile( const CurrentFile& ); // Disable copy constructor
		void operator=( const CurrentFile& ); // Disable assignment operator
	
		CurrentFile() {
            defaultPath=""; blastPath=""; testFilePath = "";
            inputDir = ""; outputDir= "";
            accnosfile = "";
            filefile = "";
            phylipfile = "";
            columnfile = "";
            listfile = "";
            rabundfile = "";
            sabundfile = "";
            namefile = "";
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
            unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
            if (concurentThreadsSupported < 1) { concurentThreadsSupported = 1; } //in case thread errors
            processors = toString(concurentThreadsSupported);
            flowfile = "";
            biomfile = "";
            countfile = "";
            summaryfile = "";
            contigsreportfile = "";
            groupMode = "group";
            sharedHeaderMode = "otu";
		}
		~CurrentFile() { instance = 0; }
};
/***********************************************/

#endif

