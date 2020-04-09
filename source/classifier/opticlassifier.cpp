//
//  opticlassifier.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/12/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "opticlassifier.hpp"

/**************************************************************************************************/
OptiClassifier::OptiClassifier(string reffasta, string reftaxonomy, string version) : Classify(){
    try {
        Utils util;
        
        //look for shortcut files - reference list file, reference otu_dist file,
        string fastaFileRoot = util.getRootName(util.getSimpleName(reffasta));
        string taxonomyFileRoot = util.getRootName(util.getSimpleName(reftaxonomy));
        
        string phyloTreeName = taxonomyFileRoot + ".tree.train";
        string refOtuDist = fastaFileRoot + ".otu.taxdist";
        
        ofstream out; ofstream out2;
        
        vector<ifstream*> files;
        ifstream* phyloTreeTest = new ifstream(phyloTreeName.c_str()); files.push_back(phyloTreeTest);
        ifstream* refDistTest = new ifstream(refOtuDist.c_str()); files.push_back(refDistTest);
 
        long start = time(NULL);
        vector< vector<string> > referenceBins;
        
        
        //if they are there make sure they were created after this release date
        bool FilesGood = false;
        if(*phyloTreeTest && *refDistTest){ FilesGood = checkReleaseDate(files, version); }

        //read taxonomy file or shortcut files to create reference otus
        if(*phyloTreeTest && *refDistTest && FilesGood){
            
            m->mothurOut("Reading template taxonomy...     "); cout.flush();
            
            phyloTree = new PhyloTree(*phyloTreeTest, phyloTreeName); //.tree.train
            maxLevel = phyloTree->getMaxLevel();
            
            m->mothurOut("DONE.\n");
            
            //m->mothurOut("Reading template distances...     "); cout.flush();
            //matrix = readOtuDistanceFile(refDistTest, referenceBins); //fills refBins
            
        }else{
            readTaxonomy(taxFile);
            
            //sanity check
            bool okay = phyloTree->ErrorCheck(names);
            
            //prevents errors caused by creating shortcut files if you had an error in the sanity check.
            if (m->getControl_pressed()) {  util.mothurRemove(phyloTreeName);  util.mothurRemove(refOtuDist); }
            else{
                m->mothurOut("Calculating template taxonomy tree...     "); cout.flush();
                
                phyloTree->printTreeNodes(phyloTreeName);
                            
                m->mothurOut("DONE.\n");
                
                //read in new phylotree with less info. - its faster
                ifstream phyloTreeTest(phyloTreeName.c_str());
                delete phyloTree;
                
                phyloTree = new PhyloTree(phyloTreeTest, phyloTreeName);
                maxLevel = phyloTree->getMaxLevel();
                
                //get reference sequences bins
                referenceBins = binReferences();
                
                //find reference otus distance
                //matrix = findReferenceOTUDistances(referenceBins);
                
            }
        }
        
        for (int i = 0; i < files.size(); i++) { delete files[i]; }

        m->mothurOut("DONE.\nIt took " + toString(time(NULL) - start) + " seconds get probabilities.\n");
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "OptiClassifier");
        exit(1);
    }
}
/**************************************************************************************************/
vector< vector<string> > OptiClassifier::binReferences() {
    try {
        
        vector< vector<string> > bins;
        vector<int> genusNodes = phyloTree->getGenusNodes();
        
        //go through nodes and build listvector
        for (int i = 0; i < genusNodes.size(); i++) {
            
             if (m->getControl_pressed()) { break; }
            
            //get parents
            TaxNode node = phyloTree->get(genusNodes[i]);
            
            vector<string> names = node.accessions;
            
            bins.push_back(names);
        }
        
        return bins;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "binReferences");
        exit(1);
    }
}

/**************************************************************************************************/

