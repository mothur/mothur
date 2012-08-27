//
//  classifysharedcommand.h
//  Mothur
//
//  Created by Abu Zaher Md. Faridee on 8/13/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__classifysharedcommand__
#define __Mothur__classifysharedcommand__

#include "command.hpp"
#include "inputdata.h"

class ClassifySharedCommand : public Command {
public:
  ClassifySharedCommand();
  ClassifySharedCommand(string option);
  ~ClassifySharedCommand() {};
  
  vector<string> setParameters();
  string getCommandName()			{ return "classify.shared";     }
  // TODO: is General is the best category for this command?
  //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
  string getCommandCategory()		{ return "General";		}
  string getOutputFileNameTag(string, string);
  string getHelpString();
  string getCitation() { return "http://www.mothur.org/wiki/classsify.shared\n"; }
  // TODO: find a proper description
  string getDescription()		{ return "find the most important otu from a shared file"; }

  int execute();
  void processSharedAndDesignData(vector<SharedRAbundVector*> lookup);
  void help() { m->mothurOut(getHelpString()); }

private:
  bool abort;
  string outputDir;
  vector<string> outputNames;
  
  string sharedfile, designfile;
  set<string> labels;
  bool allLines;
  
  GroupMap* designMap;
};

#endif /* defined(__Mothur__classifysharedcommand__) */
