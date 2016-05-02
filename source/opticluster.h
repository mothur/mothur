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

/***********************************************************************/

class OptiCluster : public Cluster {
public:
    OptiCluster(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, string, float);
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; }
    string getTag() { return("opti"); }
    
private:
    
};

#endif /* defined(__Mothur__opticluster__) */
