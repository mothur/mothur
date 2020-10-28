//
//  biomhdf5.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "biomhdf5.hpp"

/**************************************************************************************************/
BiomHDF5::BiomHDF5(string fname) : Biom(){
    try {
       
        version = "2.1.0";
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "BiomHDF5");
        exit(1);
    }
}
/**************************************************************************************************/
