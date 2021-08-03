#ifndef SPLITGROUPSCOMMAND_H
#define SPLITGROUPSCOMMAND_H

/*
 *  splitgroupscommand.h
 *  Mothur
 *
 *  Created by westcott on 9/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


/* split.groups - given a group file, split sequences and names files in to separate files *.group1.fasta and .group1.names. */


#include "command.hpp"
#include "groupmap.h"
#include "sequence.hpp"
#include "fastqread.h"
#include "getseqscommand.h"

//**********************************************************************************************************************
struct flowOutput {
    string output;
    string filename;
    int total;
    
    flowOutput(string f) { filename = f; output = ""; total = 0;  }
    flowOutput() { filename = ""; output = ""; total = 0;  }
    flowOutput(string f, string o, int t) : filename(f), output(o), total(t) {}
    
};
//**********************************************************************************************************************
struct fastqOutput {
    vector<FastqRead> output;
    string filename;
    int total;
    
    fastqOutput(string f) { filename = f;  total = 0;  }
    fastqOutput() { filename = "";  total = 0;  }
};
//**********************************************************************************************************************
struct splitGroupsStruct {
    string groupfile, countfile, namefile, inputFileName, format;
    int start, end;
    bool isFastq;
    vector<string> Groups;
    map<string, fastqOutput> parsedFastqData;
    map<string, flowOutput> parsedFlowData;
    map<string, vector<string> > outputTypes;
    vector<string> outputNames;
    MothurOut* m;
    Utils util;
    
    splitGroupsStruct(string group, string count, string name, vector<string> g, int st, int en) : groupfile(group), countfile(count), namefile(name), start(st), end(en) {
        m = MothurOut::getInstance(); format = "illumina1.8+"; isFastq = true;
        for (int i = st; i < en; i++) { Groups.push_back(g[i]); }
    }
    
    void setFiles(string input, string outfileRoot, string ext) {
        inputFileName = input;
        
        int numFlows = 0;
        if (ext != ".fastq") {
            ifstream in; util.openInputFile(input, in);
            in >> numFlows; in.close();
        }
        
        for (int i = 0; i < Groups.size(); i++) {
            string newFileName = outfileRoot + Groups[i] + ext;
            
            if (ext == ".fastq") {
                fastqOutput thisGroupsInfo(newFileName);
                parsedFastqData[Groups[i]] = thisGroupsInfo;
                ofstream out; util.openOutputFile(newFileName, out);  out.close(); //clear file for append
            }else {
                flowOutput thisGroupsInfo(newFileName);
                parsedFlowData[Groups[i]] = thisGroupsInfo;
                isFastq = false;
                
                ofstream out; util.openOutputFile(newFileName, out);  out << numFlows << endl; out.close(); //clear file for append
            }

            if (m->getControl_pressed()) { break; }
        }
    }
    
    void setFormat(string form) {  format = form;  }
};
//**********************************************************************************************************************
struct splitGroups2Struct {
    string groupfile, countfile, namefile, fastafile, listfile, outputDir;
    int start, end;
    vector<string> Groups;
    map<string, vector<string> > group2Files; //GroupName -> files(fasta, list, count) or  GroupName -> files(fasta, list, group, name) 
    map<string, vector<string> > outputTypes;
    vector<string> outputNames;
    MothurOut* m;
    Utils util;
    
    splitGroups2Struct(string group, string count, string name, vector<string> g, int st, int en) : groupfile(group), countfile(count), namefile(name), start(st), end(en) {
        m = MothurOut::getInstance();
        for (int i = st; i < en; i++) { Groups.push_back(g[i]); }
    }
    
    void setFiles(string fasta, string list, string outd) {
        fastafile = fasta;
        listfile = list;
        outputDir = outd;
        
        string fastaFileRoot = outputDir + util.getRootName(util.getSimpleName(fastafile));
        string listFileRoot = outputDir + util.getRootName(util.getSimpleName(listfile));
        string listExt = util.getExtension(listfile);
        string fastaExt = util.getExtension(fastafile);
        
        string countFileRoot, countExt, nameFileRoot, nameExt, groupFileRoot, groupExt;
        if (countfile != "") {
            countFileRoot = outputDir + util.getRootName(util.getSimpleName(countfile));
            countExt = util.getExtension(countfile);
            groupFileRoot = ""; groupExt = ""; nameFileRoot = ""; nameExt = "";
        }else {
            groupFileRoot = outputDir + util.getRootName(util.getSimpleName(groupfile));
            groupExt = util.getExtension(groupfile);
            if (namefile != "") {
                nameFileRoot = outputDir + util.getRootName(util.getSimpleName(namefile));
                nameExt = util.getExtension(namefile);
            }else { nameFileRoot = ""; nameExt = ""; }
            countFileRoot = ""; countExt = "";
        }
        
        for (int i = 0; i < Groups.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            string newListFileName = listFileRoot  + Groups[i] + listExt;
            string newFastaFileName = fastaFileRoot  + Groups[i] + fastaExt;
            string newCountFileName = countFileRoot  + Groups[i] + countExt;
            string newGroupFileName = groupFileRoot  + Groups[i] + groupExt;
            string newNameFileName = nameFileRoot  + Groups[i] + nameExt;
            
            vector<string> files;
            if (fastafile != "")    { files.push_back(newFastaFileName);    }else { files.push_back(""); }
            if (listfile != "")     { files.push_back(newListFileName);     }else { files.push_back(""); }
            
            if (countfile != "")    {
                files.push_back(newCountFileName);
            }else if (groupfile != "") {
                files.push_back(newGroupFileName);
                if (namefile != "") {
                    files.push_back(newNameFileName);
                }//else{ files.push_back(""); }
            }
            
            group2Files[Groups[i]] = files;
        }
    }
};
/***************************************************************************************/

class SplitGroupCommand : public Command {
	
public:
	SplitGroupCommand(string);	
	~SplitGroupCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "split.groups";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Split.group"; }
	string getDescription()		{ return "split a name or fasta file by group"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	vector<string> outputNames;
	vector<linePair> lines;
	string  namefile, groupfile, countfile, groups, fastafile, flowfile, fastqfile, format, listfile;
	vector<string> Groups;
	bool abort;
    int processors;
    
    void splitCountOrGroup(bool);
    void splitFastqOrFlow(string, string);
};

/***************************************************************************************/

#endif



