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

class TreeReader {
    
public:
    
    TreeReader(string tf);
	TreeReader(string tf, string gf);
    TreeReader(string tf, string gf, string nf);
	~TreeReader() {}	
    
    vector<Tree*> getTrees()            { return trees;     }
    map<string, string> getNames()      { return nameMap;   } //dups -> unique
    map<string, string> getNameMap()    { return names;     } //unique -> dups list
    
    
private:
    MothurOut* m;
	vector<Tree*> trees;
    TreeMap* tmap;
    map<string, string> nameMap; //dupName -> uniqueName
    map<string, string> names;
    
    string treefile, groupfile, namefile;
    
    bool readTrees();
    int readNamesFile();
};



#endif
