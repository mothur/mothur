#ifndef Mothur_treereader_h
#define Mothur_treereader_h

//
//  treereader.h
//  Mothur
//
//  Created by Sarah Westcott on 4/11/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "mothurout.h"
#include "tree.h"
#include "counttable.h"

class TreeReader {
    
public:
    
    TreeReader(string tf, string cf);
    TreeReader(string tf, string gf, string nf);
	~TreeReader() {}	
    
    vector<Tree*> getTrees()            { return trees;     }
    
private:
    MothurOut* m;
	vector<Tree*> trees;
    CountTable* ct;
    //map<string, string> nameMap; //dupName -> uniqueName
   // map<string, string> names;
    
    string treefile, groupfile, namefile, countfile;
    
    bool readTrees();
    int readNamesFile();
};



#endif
