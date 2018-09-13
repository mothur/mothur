//
//  optidata.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "optidata.hpp"

/***********************************************************************/
long long OptiData::print(ostream& out) {
    try {
        long long count = 0;
        for (long long i = 0; i < closeness.size(); i++) {
            out << i << '\t' << getName(i) << '\t';
            for(set<long long>::iterator it = closeness[i].begin(); it != closeness[i].end(); it++){
                out << *it << '\t';
                count++;
            }
            out << endl;
        }
        out << endl;
        return count;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "print");
        exit(1);
    }
}
/***********************************************************************/
long long OptiData::getNumClose(long long index) {
    try {
        if (index < 0) { return 0; }
        else if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return 0; }
        else { return closeness[index].size(); }
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getNumClose");
        exit(1);
    }
}
/***********************************************************************/
bool OptiData::isClose(long long i, long long toFind){
    try {
        if (i < 0) { return false; }
        else if (i > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return false; }
        
        bool found = false;
        if (closeness[i].count(toFind) != 0) { found = true; }
        return found;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "isClose");
        exit(1);
    }
}
/***********************************************************************/
set<long long> OptiData::getCloseSeqs(long long i){
    try {
        if (i < 0) { set<long long> temp; return temp; }
        else if (i > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); set<long long> temp; return temp; }
        
        return closeness[i];
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getNumClose");
        exit(1);
    }
}
/***********************************************************************/
//maps unique name to index in distance matrix
//used by sensspec to get translate the list file name to the index name for closeness shirt
map<string, long long> OptiData::getNameIndexMap() {
    try {
        map<string, long long> nameIndexes;
        for (int i = 0; i < nameMap.size(); i++) {
            vector<string> thisBinsSeqs; util.splitAtComma(nameMap[i], thisBinsSeqs);
            if (i < closeness.size()) { nameIndexes[thisBinsSeqs[0]] = i;  }
        }
        return nameIndexes;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getNameIndexMap");
        exit(1);
    }
}
/***********************************************************************/
string OptiData::getName(long long index) {
    try {
        if (index < 0) { return ""; }
        else if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return ""; }
        
        return nameMap[index];
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getName");
        exit(1);
    }
}
/***********************************************************************/
long long OptiData::getNumDists(){
    try {
        long long foundDists = 0;
        for (int i = 0; i < closeness.size(); i++) { foundDists += closeness[i].size(); }
        return foundDists;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getNumDists");
        exit(1);
    }
}
/***********************************************************************/
ListVector* OptiData::getListSingle() {
    try {
        ListVector* singlelist = NULL;
        
        if (singletons.size() == 0) { }
        else {
            singlelist = new ListVector();
            
            for (int i = 0; i < singletons.size(); i++) { singlelist->push_back(singletons[i]); }
        }
        
        return singlelist;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getListSingle");
        exit(1);
    }
}
/***********************************************************************/


