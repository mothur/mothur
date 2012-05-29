//
//  treenode.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/29/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_treenode_hpp
#define rrf_fs_prototype_treenode_hpp

  // TreeNode is a passive class, we run operations on this, but it sould not
  // contain any specific methods of "Regularized Random Forest" of itself
class TreeNode{
public:
  TreeNode(vector<TrainingSet> samples) : 
    samples(samples),
    isLeaf(false){
  }
  
  vector<TrainingSet> samples;
  bool isLeaf;
};

#endif
