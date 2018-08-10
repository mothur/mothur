//
//  fakeoptimatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/20/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "fakeoptimatrix.hpp"

/***********************************************************************/
FakeOptiMatrix::FakeOptiMatrix() : OptiData(0.03) {
    try {
        m = MothurOut::getInstance();
        
        //create 10 singletons
        for (int i = 90; i < 100; i++) { singletons.push_back(toString(i));  }
        
        //create 90 non singletons
        for (int i = 0; i < 90; i++) { nameMap.push_back(toString(i));  }
        
        closeness.resize(90);
        int count = 0;
        for (int i = 0; i < 9; i++) {
            set<int> close;
            //create list of all sequences in this set
            for (int j = 0; j < 10; j++) { close.insert((j+count)); }
            
            for (set<int>::iterator it = close.begin(); it != close.end(); it++) {
                //add close sequences to each sequence in this set, do not include self
                for (int j = 0; j < 10; j++) { if ((j+count) != *it) {   closeness[j+count].insert(*it);  } }
            }
            count += 10;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "FakeOptiMatrix", "FakeOptiMatrix");
        exit(1);
    }
}
/***********************************************************************/

