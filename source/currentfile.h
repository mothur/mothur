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
		string getPhylipFile()		{ lock_guard<std::mutex> guard(currentProtector); return phylipfile;		}
		string getColumnFile()		{ lock_guard<std::mutex> guard(currentProtector); return columnfile;		}
		string getListFile()		{ lock_guard<std::mutex> guard(currentProtector); return listfile;			}
		string getRabundFile()		{ lock_guard<std::mutex> guard(currentProtector); return rabundfile;		}
		string getSabundFile()		{ lock_guard<std::mutex> guard(currentProtector); return sabundfile;		}
		string getNameFile()		{ lock_guard<std::mutex> guard(currentProtector); return namefile;			}
		string getGroupFile()		{ lock_guard<std::mutex> guard(currentProtector); return groupfile;			}
		string getOrderFile()		{ lock_guard<std::mutex> guard(currentProtector); return orderfile;			}
		string getOrderGroupFile()	{ lock_guard<std::mutex> guard(currentProtector); return ordergroupfile;	}
		string getTreeFile()		{ lock_guard<std::mutex> guard(currentProtector); return treefile;			}
		string getSharedFile()		{ lock_guard<std::mutex> guard(currentProtector); return sharedfile;		}
		string getRelAbundFile()	{ lock_guard<std::mutex> guard(currentProtector); return relabundfile;		}
        string getCLRFile()         { lock_guard<std::mutex> guard(currentProtector); return clrfile;           }
		string getDesignFile()		{ lock_guard<std::mutex> guard(currentProtector); return designfile;		}
		string getFastaFile()		{ lock_guard<std::mutex> guard(currentProtector); return fastafile;			}
		string getSFFFile()			{ lock_guard<std::mutex> guard(currentProtector); return sfffile;			}
		string getQualFile()		{ lock_guard<std::mutex> guard(currentProtector); return qualfile;			}
		string getOligosFile()		{ lock_guard<std::mutex> guard(currentProtector); return oligosfile;		}
        string getSampleFile()      { lock_guard<std::mutex> guard(currentProtector); return samplefile;        }
        string getAccnosFile()        { lock_guard<std::mutex> guard(currentProtector); return accnosfile;        }
        string getTaxonomyFile()    { lock_guard<std::mutex> guard(currentProtector); return taxonomyfile;        }
        string getFlowFile()        { lock_guard<std::mutex> guard(currentProtector); return flowfile;            }
        string getContigsReportFile(){ lock_guard<std::mutex> guard(currentProtector); return contigsreportfile;            }
        string getBiomFile()        { lock_guard<std::mutex> guard(currentProtector); return biomfile;            }
        string getCountFile()       { lock_guard<std::mutex> guard(currentProtector); return countfile;         }
        string getSummaryFile()     { lock_guard<std::mutex> guard(currentProtector); return summaryfile;       }
        string getFileFile()        { lock_guard<std::mutex> guard(currentProtector); return filefile;          }
        string getConsTaxonomyFile(){ lock_guard<std::mutex> guard(currentProtector); return constaxonomyfile;  }
		
		void setListFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); listfile = util.getFullPathName(f);			}
        void setBiomFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); biomfile = util.getFullPathName(f);			}
        void setFlowFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); flowfile = util.getFullPathName(f);			}
        void setContigsReportFile(string f)	{ lock_guard<std::mutex> guard(currentProtector); contigsreportfile = util.getFullPathName(f);			}
        void setSummaryFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); summaryfile = util.getFullPathName(f);		}
		void setTreeFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); treefile = util.getFullPathName(f);			}
        void setGroupFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); groupfile = util.getFullPathName(f);	setGroupMode("group");		}
        void setCountFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); countfile = util.getFullPathName(f);	setGroupMode("count");		}
		void setPhylipFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); phylipfile = util.getFullPathName(f);			}
		void setColumnFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); columnfile = util.getFullPathName(f);			}
		void setNameFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); namefile = util.getFullPathName(f);			}
		void setRabundFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); rabundfile = util.getFullPathName(f);			}
		void setSabundFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); sabundfile = util.getFullPathName(f);			}
		void setSharedFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); sharedfile = util.getFullPathName(f);			}
		void setRelAbundFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); relabundfile = util.getFullPathName(f);		}
        void setCLRFile(string f)           { lock_guard<std::mutex> guard(currentProtector); clrfile = util.getFullPathName(f);        }
		void setOrderFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); orderfile = util.getFullPathName(f);			}
		void setOrderGroupFile(string f)	{ lock_guard<std::mutex> guard(currentProtector); ordergroupfile = util.getFullPathName(f);		}
		void setDesignFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); designfile = util.getFullPathName(f);			}
		void setFastaFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); fastafile = util.getFullPathName(f);			}
		void setSFFFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); sfffile = util.getFullPathName(f);			}
		void setQualFile(string f)			{ lock_guard<std::mutex> guard(currentProtector); qualfile = util.getFullPathName(f);			}
		void setOligosFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); oligosfile = util.getFullPathName(f);			}
        void setAccnosFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); accnosfile = util.getFullPathName(f);			}
        void setTaxonomyFile(string f)		{ lock_guard<std::mutex> guard(currentProtector); taxonomyfile = util.getFullPathName(f);       }
        void setConsTaxonomyFile(string f)  { lock_guard<std::mutex> guard(currentProtector); constaxonomyfile = util.getFullPathName(f);	}
        void setProgramPath(string f)       { lock_guard<std::mutex> guard(currentProtector); mothurProgramPath = util.getFullPathName(f);	}
        void setFileFile(string f)          { lock_guard<std::mutex> guard(currentProtector); filefile = util.getFullPathName(f);           }
        void setSampleFile(string f)        { lock_guard<std::mutex> guard(currentProtector); samplefile = util.getFullPathName(f);         }
    
        //current files - if you add a new type you must edit optionParser->getParameters, get.current and set.current commands and mothurOut->printCurrentFiles/clearCurrentFiles/getCurrentTypes/hasCurrentFiles. add a get and set function.
        
        string getProcessors()		{ lock_guard<std::mutex> guard(currentProtector); return processors;		}
        int setProcessors(string p);
        string getProgramPath()     { lock_guard<std::mutex> guard(currentProtector); return mothurProgramPath;  }
        string getDefaultPath() { lock_guard<std::mutex> guard(currentProtector); return defaultPath; }
        void setDefaultPath(string);
        string getTestFilePath() { lock_guard<std::mutex> guard(currentProtector); return testFilePath; }
        void setTestFilePath(string);
        string getBlastPath() { lock_guard<std::mutex> guard(currentProtector); return blastPath; }
        void setBlastPath(string);
        string getToolsPath() { lock_guard<std::mutex> guard(currentProtector); return toolsPath; }
        void setToolsPath(string);
        string getHomePath() { lock_guard<std::mutex> guard(currentProtector); return homePath; }
        void setHomePath(string);
        vector<string> getPaths() { lock_guard<std::mutex> guard(currentProtector); return paths; } //environment variable 'PATH' values
        void setPaths(vector<string>);
        string getOutputDir() { lock_guard<std::mutex> guard(currentProtector); return outputDir; }
        void setOutputDir(string f) { lock_guard<std::mutex> guard(currentProtector); outputDir = util.getFullPathName(f); }
        string getInputDir() { lock_guard<std::mutex> guard(currentProtector); return inputDir; }
        void setInputDir(string f) { lock_guard<std::mutex> guard(currentProtector); inputDir = util.getFullPathName(f); }
        void setFileName(string);
        string getReleaseDate() { lock_guard<std::mutex> guard(currentProtector); return releaseDate; }
        void setReleaseDate(string r) { lock_guard<std::mutex> guard(currentProtector); releaseDate = r; }
        string getVersion() { lock_guard<std::mutex> guard(currentProtector); return version; }
        void setVersion(string r) { lock_guard<std::mutex> guard(currentProtector); version = r; }
        vector<string> getLocations() { lock_guard<std::mutex> guard(currentProtector); vector<string> locations; locations.push_back(inputDir); locations.push_back(outputDir); locations.push_back(defaultPath); locations.push_back(mothurProgramPath); locations.push_back(toolsPath); return locations; }
    

        bool getMothurCalling()                         { lock_guard<std::mutex> guard(currentProtector); return mothurCalling;             }
        void setMothurCalling(bool t)                   { lock_guard<std::mutex> guard(currentProtector); mothurCalling = t;                }
        void printCurrentFiles(string); //string="" for just to logfile.
        void clearCurrentFiles();
        set<string> getCurrentTypes();
        bool hasCurrentFiles();
    
        string getGroupMode()                           { lock_guard<std::mutex> guard(currentProtector); return groupMode;                 }
    
        string getTestDirectory()                       { lock_guard<std::mutex> guard(currentProtector); return testDirectory;             }
        void setTestDirectory(string t)                 { lock_guard<std::mutex> guard(currentProtector); testDirectory = t;                }
   
	
	private:
		MothurOut* m;
        Utils util;
        
        vector<string> paths; //paths stored in PATH environment variables
        string logFileName, mothurProgramPath, homePath;
        string defaultPath, outputDir, blastPath, inputDir;
        string releaseDate, version;
    
        string accnosfile, phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile, taxonomyfile, biomfile, filefile, testFilePath, contigsreportfile, clrfile;
        string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile, processors, flowfile, countfile, summaryfile, constaxonomyfile, groupMode, testDirectory, sharedHeaderMode, samplefile, toolsPath;
    bool mothurCalling;
		
        void setGroupMode(string t)                     { groupMode = t;                    }
    
		static CurrentFile* instance;
		CurrentFile( const CurrentFile& ); // Disable copy constructor
		void operator=( const CurrentFile& ); // Disable assignment operator
	
        std::mutex currentProtector;
		CurrentFile() {
            m = MothurOut::getInstance();
            defaultPath=""; blastPath=""; toolsPath=""; testFilePath = "";
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
            clrfile = "";
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
            samplefile = "";
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
            mothurCalling = false;
		}
		~CurrentFile() { instance = 0; }
};
/***********************************************/

#endif

