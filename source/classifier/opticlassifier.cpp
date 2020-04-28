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
        string probFileName = fastaFileRoot + ".opti.prob";
        string probFileName2 = fastaFileRoot + ".opti.numNonZero";
        
        ofstream out; ofstream out2;
        
        vector<ifstream*> files;
        ifstream* phyloTreeTest = new ifstream(phyloTreeName.c_str()); files.push_back(phyloTreeTest);
        ifstream* probFileTest2 = new ifstream(probFileName2.c_str()); files.push_back(probFileTest2);
        ifstream* probFileTest = new ifstream(probFileName.c_str());   files.push_back(probFileTest);
        
        long start = time(NULL);
        
        //if they are there make sure they were created after this release date
        bool FilesGood = false;
        if(*phyloTreeTest && *probFileTest && *probFileTest2){ FilesGood = checkReleaseDate(files, version); }

        //read taxonomy file or shortcut files to create reference otus
        if(*phyloTreeTest && *probFileTest && *probFileTest2 && FilesGood){
            
            m->mothurOut("Reading template taxonomy...     "); cout.flush();
            
            phyloTree = new PhyloTree(*phyloTreeTest, phyloTreeName); //.tree.train
            maxLevel = phyloTree->getMaxLevel();
            
            m->mothurOut("DONE.\n");
            
            //m->mothurOut("Reading template distances...     "); cout.flush();
            //matrix = readOtuDistanceFile(refDistTest, referenceBins); //fills refBins
            
        }else{
            generateDatabaseAndNames(reftaxonomy, reffasta, "opti", 8, 0.0, 0.0, 0.0, 0.0, version);
            
            //sanity check
            bool okay = phyloTree->ErrorCheck(names);
            
            //prevents errors caused by creating shortcut files if you had an error in the sanity check.
            if (m->getControl_pressed()) {  util.mothurRemove(phyloTreeName);  util.mothurRemove(probFileName); util.mothurRemove(probFileName2); }
            else{
                genusNodes = phyloTree->getGenusNodes();
                genusTotals = phyloTree->getGenusTotals();
                
                m->mothurOut("Calculating template taxonomy tree...     "); cout.flush();
                
                phyloTree->printTreeNodes(phyloTreeName);
                            
                m->mothurOut("DONE.\n");
                
                int numAlignedColumns = database->getLongestBase();
                int numReferences = names.size();
                charGenusProb.resize(numAlignedColumns);
                map<char, int> baseMap; map<char, int>::iterator itBase;
                baseMap['A'] = 0;
                baseMap['T'] = 1;
                baseMap['G'] = 2;
                baseMap['C'] = 3;
                baseMap['-'] = 4;
                
                if (shortcuts) {
                    util.openOutputFile(probFileName, out);
                
                    //output mothur version
                    out << "#" << version << endl;
                
                    out << numAlignedColumns << endl;
                
                    util.openOutputFile(probFileName2, out2);
                
                    //output mothur version
                    out2 << "#" << version << endl;
                }
                
                //for each column in the alignment
                for (int i = 0; i < numAlignedColumns; i++) {
                    
                    if (shortcuts) {  out << i << '\t'; }
                    
                    if (m->getControl_pressed()) {  break; }
                    
                    char allSame = 'x';
                    vector< vector<int> > thisColumnDistribution = database->get(i, allSame);
                    
                    vector<int> numNotZero;  numNotZero.resize(5, 0);
                    if (allSame != 'x') { //all bases are the same in this location, simplify calc, ignore thisColumnDistribution
                        
                        //probabilityInTemplate = (# of seqs with that char at this location in template + 0.50) / (total number of seqs in template + 1);
                        
                        //probabilityInThisTaxonomy = (# of seqs with that char at this location in this taxonomy + probabilityInTemplate) / (total number of seqs in this taxonomy + 1);
                        for (itBase = baseMap.begin(); itBase != baseMap.end(); itBase++) {
                            
                            float probabilityInTemplate = 0;
                            for (int k = 0; k < genusNodes.size(); k++) {
                                
                                if (itBase->first == allSame) {
                                    probabilityInTemplate = (numReferences + 0.50) / (float) (numReferences + 1);
                                    charGenusProb[i][baseMap[allSame]][k] = log((genusTotals[k] + probabilityInTemplate) / (float) (genusTotals[k] + 1));
                                    
                                    if (shortcuts) { out << baseMap[allSame] << '\t' << k << '\t' << charGenusProb[i][baseMap[allSame]][k] << '\t' ; }
                                    
                                    numNotZero[baseMap[allSame]]++;
                                }else {
                                    probabilityInTemplate = (0 + 0.50) / (float) (numReferences + 1);
                                    charGenusProb[i][itBase->second][k] = log((0 + probabilityInTemplate) / (float) (genusTotals[k] + 1));
                                }
                            }
                            
                            if (shortcuts) {
                                out << endl;
                                out2 << itBase->second << '\t' << probabilityInTemplate << '\t' << numNotZero[itBase->second] << '\t' << endl;
                            }
                        }
                    }else {
                        
                        for (int j = 0; j < thisColumnDistribution.size(); j++) { //should be 5, A,T,G,C,-
                            //probabilityInTemplate = (# of seqs with that char at this location in template + 0.50) / (total number of seqs in template + 1);
                            float probabilityInTemplate = (thisColumnDistribution[j].size() + 0.50) / (float) (numReferences + 1);
                            
                            vector<int> count; count.resize(genusNodes.size(), 0);
                            for (int k = 0; k < thisColumnDistribution[j].size(); k++) { //find genus of seqs
                                int temp = phyloTree->getGenusIndex(names[thisColumnDistribution[j][k]]);
                                count[temp]++;  //increment count of seq in this genus who have this base at this location
                            }
                            
                            for (int k = 0; k < genusNodes.size(); k++) {
                                charGenusProb[i][j][k] = log((count[k] + probabilityInTemplate) / (float) (genusTotals[k] + 1));
                                
                                if (count[k] != 0) {
                                    if (shortcuts) { out << j << '\t' << k << '\t' << charGenusProb[i][j][k] << '\t' ; }
                                    numNotZero[j]++;
                                }
                            }
                            
                            if (shortcuts) {
                                out << endl;
                                out2 << j << '\t' << probabilityInTemplate << '\t' << numNotZero[j] << '\t' << endl;
                            }
                        }
                    }
                    
                    
                }
                if (shortcuts) { out.close(); out2.close();  }
                
                //read in new phylotree with less info. - its faster
                ifstream phyloTreeTest(phyloTreeName.c_str());
                delete phyloTree;
                
                phyloTree = new PhyloTree(phyloTreeTest, phyloTreeName);
                maxLevel = phyloTree->getMaxLevel();
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
void OptiClassifier::readProbFile(ifstream& in, ifstream& inNum) {
    try{
        Utils util;
        //read version
        string line = util.getline(in); util.gobble(in);
        
        int numColumnsAlignment = 0;
        in >> numColumnsAlignment; util.gobble(in);
    
        //initialze probabilities
        charGenusProb.resize(numColumnsAlignment);
        for (int j = 0; j < charGenusProb.size(); j++) {
            charGenusProb[j].resize(5);
            for (int i = 0; i < 5; i++) {
                charGenusProb[j][i].resize(genusNodes.size());
            }
        }
        /*
        int name, count;  count = 0;
        vector<int> num; num.resize(numColumnsAlignment);
        float prob;
        vector<float> zeroCountProb; zeroCountProb.resize(numColumnsAlignment);
        
        //read version
        string line2 = util.getline(inNum); util.gobble(inNum);
        float probTemp;
        
        while (inNum) {
            inNum >> zeroCountProb[count] >> num[count] >> probTemp;
            count++;
            util.gobble(inNum);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(zeroCountProb[count]) + '\t' + toString(num[count]) + '\t' + toString(numKmers) + "\n"); }

            
        }
        inNum.close();
 
        while(in) {
            in >> kmer;
            
            //set them all to zero value
            for (int i = 0; i < genusNodes.size(); i++) {
                wordGenusProb[kmer][i] = log(zeroCountProb[kmer] / (float) (genusTotals[i]+1));
            }
           
            //get probs for nonzero values
            for (int i = 0; i < num[kmer]; i++) {
                in >> name >> prob;
                wordGenusProb[kmer][name] = prob;
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(name) + '\t' + toString(prob) + '\t' + toString(kmer) + "\n"); }
            }
            
            util.gobble(in);
        }
        in.close();
          */
    }
    catch(exception& e) {
        m->errorOut(e, "Bayesian", "readProbFile");
        exit(1);
    }
}
/**************************************************************************************************/

