//
//  biomhdf5.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef biomhdf5_hpp
#define biomhdf5_hpp

#include "biom.hpp"

//http://biom-format.org/documentation/format_versions/biom-2.1.html

class BiomHDF5 : public Biom {
    
public:
    
    BiomHDF5(string);
    ~BiomHDF5() {  }
    
    
    
private:
   
    
};

#endif /* biomhdf5_hpp */
