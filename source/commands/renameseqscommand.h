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
	RenameSeqsCommand();
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
    
	string fastaFile, nameFile, groupfile, outputDir, placement, delim, countfile, qualfile, contigsfile, fileFile, mapFile;
	vector<string> outputNames;
	bool abort;
	
	map<string, string> nameMap;
    int readQual(map<string, string>&);
    int readContigs(map<string, string>&);
    int readFasta(string, map<string, string>&);
    int processFile(map<string, string>&);
    int readMapFile(map<string, string>&);
    vector< map<string, string> > readFiles();
    
};



#endif
