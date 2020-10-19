//
//  renameseqscommand.h
//  Mothur
//
//  Created by SarahsWork on 5/28/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_renameseqscommand_h
#define Mothur_renameseqscommand_h

#include "command.hpp"
#include "filefile.hpp"

class RenameSeqsCommand : public Command {
    
#ifdef UNIT_TEST
    friend class TestRenameSeqsCommand;
#endif
    
public:
	RenameSeqsCommand(string);
	~RenameSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "rename.seqs";		}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();
    string getOutputPattern(string);
	string getCitation() { return "http://www.mothur.org/wiki/Rename.seqs"; }
	string getDescription()		{ return "rename sequences"; }
    
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
	
private:
    
	string fastaFile, listfile, nameFile, groupfile,  placement, delim, countfile, qualfile, contigsfile, fileFile, mapFile, taxfile;
	vector<string> outputNames;
	bool abort, ignoreNew;
	
	map<string, string> nameMap;
    void readQual(map<string, string>&);
    void readTax(map<string, string>&);
    void readContigs(map<string, string>&);
    void readList(map<string, string>&);
    void readFasta(map<string, string>&);
    int processFile(map<string, string>&);
    int readMapFile(map<string, string>&);
    vector< map<string, string> > readFiles();
    void processNameGroupCountFiles(map<string, string>&, map<string, string>&);
    
};



#endif
