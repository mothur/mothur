//
//  opticluster.h
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__opticluster__
#define __Mothur__opticluster__

#include "cluster.hpp"
#include "optimatrix.h"

/***********************************************************************/

class OptiCluster : public Cluster {
    
public:
    OptiCluster(OptiMatrix*);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { return("opti"); }
    
    int initialize();  //randomize and place in "best" OTUs
    bool update(double&); //returns whether list changed and MCC
    
private:
    ListVector* list;
    
    bool updateSplit(double oldMCC, double& newMCC); //return changed
    bool updateMerge(double oldMCC, double& newMCC); //returns changed
    double calcMCC(double, double, double, double);
    
};

#endif /* defined(__Mothur__opticluster__) */
