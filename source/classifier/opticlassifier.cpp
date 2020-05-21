//
//  opticlassifier.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/12/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "opticlassifier.hpp"

/**************************************************************************************************/
OptiClassifier::OptiClassifier(string reftaxonomy, string reffasta, int cutoff, int i, bool sh, string version) : Classify(), confidenceThreshold(cutoff), iters(i) {
    try {
        Utils util;
        reverse['A'] = 'T';
        reverse['T'] = 'A';
        reverse['C'] = 'G';
        reverse['G'] = 'C';
        reverse['-'] = '-';
        reverse['N'] = 'N';
        numBases = 6;
        shortcuts = sh;
        
        string baseName = reftaxonomy;
        string baseTName = reffasta;
           
        string tfileroot = util.getFullPathName(baseTName.substr(0,baseTName.find_last_of(".")+1));
        string tempfileroot = util.getRootName(util.getSimpleName(baseName));
        string phyloTreeName = tfileroot + "tree.train";
        string probFileName = tfileroot + tempfileroot + "opti.prob";
        string probFileName2 = tfileroot + tempfileroot + "opti.numNonZero";
                        
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
            
            phyloTree = new PhyloTree(*phyloTreeTest, phyloTreeName);
            maxLevel = phyloTree->getMaxLevel();
            
            m->mothurOut("DONE.\n");
            
            genusNodes = phyloTree->getGenusNodes();
            genusTotals = phyloTree->getGenusTotals();
            
            m->mothurOut("Reading template probabilities...     "); cout.flush();
            readProbFile(*probFileTest, *probFileTest2);
            
        }else{
            generateDatabaseAndNames(reftaxonomy, reffasta, "opti", 8, 0.0, 0.0, 0.0, 0.0, version);
            
            //prevents errors caused by creating shortcut files if you had an error in the sanity check.
            if (m->getControl_pressed()) {  util.mothurRemove(phyloTreeName);  util.mothurRemove(probFileName); util.mothurRemove(probFileName2); }
            else{
                genusNodes = phyloTree->getGenusNodes();
                genusTotals = phyloTree->getGenusTotals();
                
                m->mothurOut("Calculating template taxonomy tree...     "); cout.flush();
                
                phyloTree->printTreeNodes(phyloTreeName);
                            
                m->mothurOut("DONE.\n");
                
                numAlignedColumns = database->getLongestBase();
                int numReferences = names.size();
                
                baseMap['A'] = 0; mapBase[0] = 'A';
                baseMap['T'] = 1; mapBase[1] = 'T';
                baseMap['G'] = 2; mapBase[2] = 'G';
                baseMap['C'] = 3; mapBase[3] = 'C';
                baseMap['-'] = 4; mapBase[4] = '-';
                baseMap['N'] = 5; mapBase[5] = 'N';
                
                baseProbs['A'] = 0.0;
                baseProbs['T'] = 0.0;
                baseProbs['G'] = 0.0;
                baseProbs['C'] = 0.0;
                baseProbs['-'] = 0.0;
                baseProbs['N'] = 0.0;
                
                //initailize charGenusProb[alignmentLength][numBases][numberOfGenus]
                charGenusProb.resize(numAlignedColumns);
                reversedProbs.resize(numAlignedColumns, baseProbs);
                
                for (int i = 0; i < numAlignedColumns; i++) {
                    charGenusProb[i].resize(numBases);
                    for (int j = 0; j < numBases; j++) { charGenusProb[i][j].resize(genusNodes.size(), 0.0); }
                }

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
                    
                    int numNotZero = 0;
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
                                    
                                    numNotZero++;
                                }else { //zero count prob
                                    probabilityInTemplate = (0 + 0.50) / (float) (numReferences + 1);
                                    charGenusProb[i][itBase->second][k] = log((0 + probabilityInTemplate) / (float) (genusTotals[k] + 1));
                                }
                            }
                            
                            //save to determine if sequence is reversed
                            reversedProbs[i][itBase->second] = probabilityInTemplate;
                            
                            if (shortcuts) {
                                out2 << itBase->second << '\t' << probabilityInTemplate << '\t';
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
                                    numNotZero++;
                                }
                            }
                            
                            //save to determine if sequence is reversed
                            reversedProbs[i][mapBase[j]] = probabilityInTemplate;
                            
                            if (shortcuts) {
                                out2 << j << '\t' << probabilityInTemplate << '\t';
                            }
                        }
                    }
                    if (shortcuts) {
                        out << endl;
                        out2 << numNotZero << endl;
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

        for (int i = 0; i < numAlignedColumns; i++) { allCols.push_back(i); }
        
        m->mothurOut("DONE.\nIt took " + toString(time(NULL) - start) + " seconds get probabilities.\n");
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "OptiClassifier");
        exit(1);
    }
}
/**************************************************************************************************/
string OptiClassifier::getTaxonomy(Sequence* seq, string& simpleTax, bool& flipped) {
    try {
        string tax = "";
        simpleTax = "";
        flipped = false;
        
        string aligned = seq->getAligned();
        
        if (aligned.length() != numAlignedColumns) {  m->mothurOut("[ERROR]: Alignment mismatch, cannot classify "  + seq->getName() + ". The reference sequences have an aligned length of " + toString(numAlignedColumns) + " but the length of this sequence is " + toString(aligned.length()) + ".\n");  simpleTax = "unknown;";  return "unknown;"; }
        
        //convert '.' gaps to '-'
        for (int i = 0; i < numAlignedColumns; i++) { if (aligned[i] == '.') { aligned[i] = '-'; } }
        
        if (isReversed(aligned)) {
            flipped = true;
            seq->reverseComplement();
            aligned = seq->getAligned();
        }

        int index = getMostProbableTaxonomy(aligned, allCols);
        
        TaxNode taxonomyTemp = phyloTree->get(index);
        
        string thisTax = "";
        while (taxonomyTemp.level != 0) { //while you are not at the root
            thisTax = taxonomyTemp.name + ";" + thisTax;
            taxonomyTemp = phyloTree->get(taxonomyTemp.parent);
        }
        if (m->getControl_pressed()) { return tax; }
    
        if (m->getDebug()) {  m->mothurOut(seq->getName() + "\t"); }
        
        //bootstrap - to set confidenceScore
        tax = bootstrapResults(aligned, index, simpleTax);
        
        if (m->getDebug()) {  m->mothurOut("\n"); }
        
        return tax;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "getTaxonomy");
        exit(1);
    }
}
/**************************************************************************************************/
string OptiClassifier::bootstrapResults(string aligned, int tax, string& simpleTax) {
    try {
                
        map<int, int> confidenceScores;
        
        //initialize confidences to 0
        int seqIndex = tax;
        TaxNode seq = phyloTree->get(tax);
        confidenceScores[tax] = 0;
        
        while (seq.level != 0) { //while you are not at the root
            seqIndex = seq.parent;
            confidenceScores[seqIndex] = 0;
            seq = phyloTree->get(seq.parent);
        }
                
        map<int, int>::iterator itBoot;
        map<int, int>::iterator itBoot2;
        map<int, int>::iterator itConvert;
        
        Utils util;
        for (int i = 0; i < iters; i++) {
            if (m->getControl_pressed()) { return "control"; }
            
            vector<int> sampledCols;
            for (int j = 0; j < numAlignedColumns; j++) {
                int index = util.getRandomIndex(numAlignedColumns-1);
            
                sampledCols.push_back(index);
            }
            
            //get taxonomy
            int newTax = getMostProbableTaxonomy(aligned, sampledCols);
            //int newTax = 1;
            TaxNode taxonomyTemp = phyloTree->get(newTax);
            string thisTax = "";

            //add to confidence results
            while (taxonomyTemp.level != 0) { //while you are not at the root
                itBoot2 = confidenceScores.find(newTax); //is this a classification we already have a count on
                
                if (itBoot2 != confidenceScores.end()) { //this is a classification we need a confidence for
                    (itBoot2->second)++;
                }
                
                thisTax = taxonomyTemp.name + ";" + thisTax;
                newTax = taxonomyTemp.parent;
                taxonomyTemp = phyloTree->get(newTax);
            }
        }
        
        string confidenceTax = "";
        simpleTax = "";
        
        int seqTaxIndex = tax;
        TaxNode seqTax = phyloTree->get(tax);
        
        while (seqTax.level != 0) { //while you are not at the root
                    
                itBoot2 = confidenceScores.find(seqTaxIndex); //is this a classification we already have a count on
                
                int confidence = 0;
                if (itBoot2 != confidenceScores.end()) { //already in confidence scores
                    confidence = itBoot2->second;
                }
                
                if (m->getDebug()) { m->mothurOut(seqTax.name + "(" + toString(((confidence/(float)iters) * 100)) + ");"); }
            
                if (((confidence/(float)iters) * 100) >= confidenceThreshold) {
                    confidenceTax = seqTax.name + "(" + toString(((confidence/(float)iters) * 100)) + ");" + confidenceTax;
                    simpleTax = seqTax.name + ";" + simpleTax;
                }
                
                
                seqTaxIndex = seqTax.parent;
                seqTax = phyloTree->get(seqTax.parent);
        }
        
        if (confidenceTax == "") { confidenceTax = "unknown;"; simpleTax = "unknown;";  }
    
        return confidenceTax;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "bootstrapResults");
        exit(1);
    }
}
/**************************************************************************************************/
int OptiClassifier::getMostProbableTaxonomy(string aligned, vector<int> cols) {
    try {
        int indexofGenus = 0;
        
        double maxProbability = -1000000.0;
        //find taxonomy with highest probability that this sequence is from it

        for (int k = 0; k < genusNodes.size(); k++) {
            //for each taxonomy calc its probability
            
            double prob = 0.0000;
            for (int i = 0; i < cols.size(); i++) {
                int indexOfBaseInColsI = baseMap[aligned[cols[i]]]; //0 -> A, 1 -> T, 2 -> C, 3 -> G, 4 -> -, 5 -> N
                prob += charGenusProb[cols[i]][indexOfBaseInColsI][k];
            }

            //is this the taxonomy with the greatest probability?
            if (prob > maxProbability) {
                indexofGenus = genusNodes[k];
                maxProbability = prob;
            }
        }
            
        return indexofGenus;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "getMostProbableTaxonomy");
        exit(1);
    }
}
//********************************************************************************************************************
//if it is more probable that the reverse compliments are in the template, then we assume the sequence is reversed.
//vector< map< char, float> > reversedProbs; //reversedProbs[0]['A'] = Probability of 'A' being at alignment location 0 in the reference. reversedProbs[0]['T'] = Probability of 'T' being at alignment location 0 in the reference.
bool OptiClassifier::isReversed(string aligned){
    try{
        bool reversed = false;
        float prob = 0;
        float reverseProb = 0;
         
        for (int i = 0; i < numAlignedColumns; i++){
            int base = aligned[i];
            
            prob += reversedProbs[i][base];
            reverseProb += reversedProbs[i][reverse[base]];
        }
        
        if (reverseProb > prob){ reversed = true; }
    
        return reversed;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "isReversed");
        exit(1);
    }
}
/**************************************************************************************************/
void OptiClassifier::readProbFile(ifstream& in, ifstream& inNum) {
    try{
        Utils util;
        //read version
        string line = util.getline(in); util.gobble(in);
        
        in >> numAlignedColumns; util.gobble(in);
    
        //initialze probabilities
        int numGenus = genusNodes.size();
        charGenusProb.resize(numAlignedColumns);
        for (int i = 0; i < charGenusProb.size(); i++) {
            charGenusProb[i].resize(numBases);
            for (int j = 0; j < numBases; j++) {
                charGenusProb[i][j].resize(numGenus);
            }
        }
        
        baseProbs['A'] = 0.0;
        baseProbs['T'] = 0.0;
        baseProbs['G'] = 0.0;
        baseProbs['C'] = 0.0;
        baseProbs['-'] = 0.0;
        baseProbs['N'] = 0.0;
        
        baseMap['A'] = 0; mapBase[0] = 'A';
        baseMap['T'] = 1; mapBase[1] = 'T';
        baseMap['G'] = 2; mapBase[2] = 'G';
        baseMap['C'] = 3; mapBase[3] = 'C';
        baseMap['-'] = 4; mapBase[4] = '-';
        baseMap['N'] = 5; mapBase[5] = 'N';
        
        reversedProbs.resize(numAlignedColumns, baseProbs);
        
        int base, alignmentLocation;  alignmentLocation = 0;
        vector<int> num; num.resize(numAlignedColumns); //num nonzero probs for this alignment location
        vector< vector<float> > probabilityInTemplate; probabilityInTemplate.resize(numAlignedColumns);
        for (int i = 0; i < numAlignedColumns; i++) { probabilityInTemplate[i].resize(numBases, 0); }
        
        //read version
        string line2 = util.getline(inNum); util.gobble(inNum);
        
        while (!inNum.eof()) {
            
            for (int i = 0; i < numBases; i++) {
                inNum >> base; util.gobble(inNum);
                inNum >> probabilityInTemplate[alignmentLocation][base]; util.gobble(inNum);
                
                reversedProbs[alignmentLocation][mapBase[base]] = probabilityInTemplate[alignmentLocation][base];
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(base) + '\t' + toString(probabilityInTemplate[alignmentLocation][base]) + '\t'); }
            }
            inNum >> num[alignmentLocation]; util.gobble(inNum);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(num[alignmentLocation])  + "\n"); }
            
            alignmentLocation++;
        }
        inNum.close();

        while(!in.eof()) {
            in >> alignmentLocation;
            
            for (int j = 0; j < numBases; j++) {
                
                //set them all to zero value
                for (int k = 0; k < genusNodes.size(); k++) {
                    charGenusProb[alignmentLocation][j][k] = log(probabilityInTemplate[alignmentLocation][j] / (float) (genusTotals[k]+1));
                }
            }
            
            //update non zero values
            int genus = 0; double prob = 0.0;
            for (int i = 0; i < num[alignmentLocation]; i++) {
                in >> base >> genus >> prob;
                charGenusProb[alignmentLocation][base][genus] = prob;
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(alignmentLocation) + '\t' + toString(base) + '\t' + toString(genus) + '\t' + toString(prob) + "\n"); }
            }
             util.gobble(in);
        }
        in.close();
    }
    catch(exception& e) {
        m->errorOut(e, "OptiClassifier", "readProbFile");
        exit(1);
    }
}
/**************************************************************************************************/

