//
//  kimura.hpp
//  Mothur
//
//  Created by Sarah Westcott on 7/2/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef kimura_hpp
#define kimura_hpp

#include "calculator.h"


/**************************************************************************************************/
//Kimura(Kimura formula)

class Kimura : public DistCalc {
    
public:
    
    Kimura(double c) : DistCalc(c) { name = "Kimura (Kimura formula)"; }

    double calcDist(Protein A, Protein B); //calc distance between 2 seqeunces
    string getCitation() { return "https://evolution.gs.washington.edu/phylip/doc/protdist.html, https://evolution.genetics.washington.edu/phylip/credits.html"; }
    
private:
    
    
};

/**************************************************************************************************/

#endif /* kimura_hpp */
