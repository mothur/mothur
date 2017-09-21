//
//  kmerTree.cpp
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#include "kmernode.h"
#include "kmertree.h"

/**************************************************************************************************/

KmerTree::KmerTree(string referenceFileName, string taxonomyFileName, int k, int cutoff) : Classify(), confidenceThreshold(cutoff), kmerSize(k){
	try {
        KmerNode* newNode = new KmerNode("Root", 0, kmerSize);
        tree.push_back(newNode);			//	the tree is stored as a vector of elements of type TaxonomyNode
        
        int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
        numPossibleKmers = power4s[kmerSize];
        
        string refTaxonomy;
        
        readTaxonomy(taxonomyFileName);
        
        ifstream referenceFile;
        m->openInputFile(referenceFileName, referenceFile);
        bool error = false;
        while(!referenceFile.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(referenceFile);  m->gobble(referenceFile);
            
            if (seq.getName() != "") {
                map<string, string>::iterator it = taxonomy.find(seq.getName());
                
                if (it != taxonomy.end()) {
                    refTaxonomy = it->second;		//	lookup the taxonomy string for the current reference sequence
                    vector<int> kmerProfile = ripKmerProfile(seq.getUnaligned());	//convert to kmer vector
                    addTaxonomyToTree(seq.getName(), refTaxonomy, kmerProfile);
                }else {
                    m->mothurOut(seq.getName() + " is in your reference file, but not in your taxonomy file, please correct.\n"); error = true;
                }
            }
        }
        referenceFile.close();
        
        if (error) { m->setControl_pressed(true); }
        
        numTaxa = (int)tree.size();
        numLevels = 0;
        for(int i=0;i<numTaxa;i++){
            int level = tree[i]->getLevel();
            if(level > numLevels){	numLevels = level;	}
        }
        numLevels++;
        
        aggregateThetas();
        
        int dbSize = tree[0]->getNumSeqs();
        
        for(int i=0;i<numTaxa;i++){
            tree[i]->checkTheta();
            tree[i]->setNumUniqueKmers(tree[0]->getNumUniqueKmers());
            tree[i]->setTotalSeqs(dbSize);
        }
    }
	catch(exception& e) {
		m->errorOut(e, "KmerTree", "KmerTree");
		exit(1);
	}
}

/**************************************************************************************************/

KmerTree::~KmerTree(){
	
	for(int i=0;i<tree.size();i++){
		delete tree[i];
	}
	
}	
/**********************************************************************************************************************/

vector<int> KmerTree::ripKmerProfile(string sequence){
    try {
        //	assume all input sequences are unaligned
        
        int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
        
        int nKmers = (int)sequence.length() - kmerSize + 1;
        
        vector<int> kmerProfile(numPossibleKmers + 1, 0);
        
        for(int i=0;i<nKmers;i++){
            
            if (m->getControl_pressed()) { break; }
            
            int kmer = 0;
            for(int j=0;j<kmerSize;j++){
                if(toupper(sequence[j+i]) == 'A')		{	kmer += (0 * power4s[kmerSize-j-1]);	}
                else if(toupper(sequence[j+i]) == 'C')	{	kmer += (1 * power4s[kmerSize-j-1]);	}
                else if(toupper(sequence[j+i]) == 'G')	{	kmer += (2 * power4s[kmerSize-j-1]);	}
                else if(toupper(sequence[j+i]) == 'U')	{	kmer += (3 * power4s[kmerSize-j-1]);	}
                else if(toupper(sequence[j+i]) == 'T')	{	kmer += (3 * power4s[kmerSize-j-1]);	}
                else									{	kmer = power4s[kmerSize]; j = kmerSize;	}
            }
            kmerProfile[kmer] = 1;
        }
        
        return kmerProfile;	
    }
	catch(exception& e) {
		m->errorOut(e, "KmerTree", "ripKmerProfile");
		exit(1);
	}
}

/**************************************************************************************************/

int KmerTree::addTaxonomyToTree(string seqName, string taxonomy, vector<int>& sequence){
	try {
        KmerNode* newNode;
        string taxonName = "";
        int treePosition = 0;							//	the root is element 0
        
        
        int level = 1;
        
        for(int i=0;i<taxonomy.length();i++){			//	step through taxonomy string...
            
            if (m->getControl_pressed()) { break; }
            if(taxonomy[i] == ';'){						//	looking for semicolons...
                
                if (taxonName == "") {  m->mothurOut(seqName + " has an error in the taxonomy.  This may be due to a ;;"); m->mothurOutEndLine(); m->setControl_pressed(true); }
                
                int newIndex = tree[treePosition]->getChildIndex(taxonName);//	look to see if your current node already
                //	   has a child with the new taxonName
                if(newIndex != -1)	{	treePosition = newIndex;	}		//	if you've seen it before, jump to that
                else {														//	   position in the tree
                    int newChildIndex = (int)tree.size();					//	otherwise, we'll have to create one...
                    tree[treePosition]->makeChild(taxonName, newChildIndex);
                    
                    newNode = new KmerNode(taxonName, level, kmerSize);
                    
                    newNode->setParent(treePosition);
                    
                    tree.push_back(newNode);
                    treePosition = newChildIndex;
                }
                
                //	sequence data to that node to update that node's theta - seems slow...				
                taxonName = "";								//	clear out the taxon name that we will build as we look 
                level++;
                
            }												//	for a semicolon
            else{
                taxonName += taxonomy[i];					//	keep adding letters until we reach a semicolon
            }
        }
        
        tree[treePosition]->loadSequence(sequence);	//	now that we've gotten to the correct node, add the
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "KmerTree", "addTaxonomyToTree");
		exit(1);
	}
	
}

/**************************************************************************************************/

int KmerTree::aggregateThetas(){
	try {
        vector<vector<int> > levelMatrix(numLevels+1);
        
        for(int i=0;i<tree.size();i++){
            if (m->getControl_pressed()) { return 0; }
            levelMatrix[tree[i]->getLevel()].push_back(i);
        }
        
        for(int i=numLevels-1;i>0;i--) {
            if (m->getControl_pressed()) { return 0; }
            
            for(int j=0;j<levelMatrix[i].size();j++){
                
                KmerNode* holder = tree[levelMatrix[i][j]];
                
                tree[holder->getParent()]->addThetas(holder->getTheta(), holder->getNumSeqs());				
            }
        }
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "KmerTree", "aggregateThetas");
		exit(1);
	}
}

/**************************************************************************************************/

int KmerTree::getMinRiskIndexKmer(vector<int>& sequence, vector<int>& taxaIndices, vector<double>& probabilities){
	try {
        int numProbs = (int)probabilities.size();
        
        vector<double> G(numProbs, 0.2);	//a random sequence will, on average, be 20% similar to any other sequence; not sure that this holds up for kmers; whatever.
        vector<double> risk(numProbs, 0);
        
        for(int i=1;i<numProbs;i++){ //use if you want the outlier group
            if (m->getControl_pressed()) { return 0; }
            G[i] = tree[taxaIndices[i]]->getSimToConsensus(sequence);
        }
        
        double minRisk = 1e6;
        int minRiskIndex = 0;
        
        for(int i=0;i<numProbs;i++){
            if (m->getControl_pressed()) { return 0; }
            for(int j=0;j<numProbs;j++){
                if(i != j){
                    risk[i] += probabilities[j] * G[j];
                }			
            }
            
            if(risk[i] < minRisk){
                minRisk = risk[i];
                minRiskIndex = i;
            }
        }
        
        return minRiskIndex;
    }
	catch(exception& e) {
		m->errorOut(e, "KmerTree", "getMinRiskIndexKmer");
		exit(1);
	}
}

/**************************************************************************************************/

int KmerTree::sanityCheck(vector<vector<int> >& indices, vector<int>& maxIndices){
	try {
        int finalLevel = (int)indices.size()-1;
        
        for(int position=1;position<indices.size();position++){
            if (m->getControl_pressed()) { return 0; }
            int predictedParent = tree[indices[position][maxIndices[position]]]->getParent();
            int actualParent = indices[position-1][maxIndices[position-1]];
            
            if(predictedParent != actualParent){
                finalLevel = position - 1;
                return finalLevel;
            }
        }
        return finalLevel;
	}
	catch(exception& e) {
		m->errorOut(e, "KmerTree", "sanityCheck");
		exit(1);
	}
}

/**************************************************************************************************/
string KmerTree::getTaxonomy(Sequence* thisSeq, string& simpleTax, bool& flipped){
	try {
        simpleTax = "";
        string seqName = thisSeq->getName(); string querySequence = thisSeq->getAligned(); string taxonProbabilityString = "";
        string unalignedSeq = thisSeq->getUnaligned();
        
        double logPOutlier = (querySequence.length() - kmerSize + 1) * log(1.0/(double)tree[0]->getNumUniqueKmers());
        
        vector<int> queryProfile = ripKmerProfile(unalignedSeq);	//convert to kmer vector
        
        vector<vector<double> > pXgivenKj_D_j(numLevels);
        vector<vector<int> > indices(numLevels);
        for(int i=0;i<numLevels;i++){
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            pXgivenKj_D_j[i].push_back(logPOutlier);
            indices[i].push_back(-1);
        }
        
        for(int i=0;i<numTaxa;i++){
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            pXgivenKj_D_j[tree[i]->getLevel()].push_back(tree[i]->getPxGivenkj_D_j(queryProfile));
            indices[tree[i]->getLevel()].push_back(i);
        }
        
        vector<double> sumLikelihood(numLevels, 0);
        vector<double> bestPosterior(numLevels, 0);
        vector<int> maxIndex(numLevels, 0);
        int maxPosteriorIndex;
        
        //let's find the best level and taxa within that level
        for(int i=0;i<numLevels;i++){ //go across all j's - from the root to genus
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            
            int numTaxaInLevel = (int)indices[i].size();
            
            vector<double> posteriors(numTaxaInLevel, 0);		
            sumLikelihood[i] = getLogExpSum(pXgivenKj_D_j[i], maxPosteriorIndex);
            
            maxPosteriorIndex = 0;
            for(int j=0;j<numTaxaInLevel;j++){
                posteriors[j] = exp(pXgivenKj_D_j[i][j] - sumLikelihood[i]);
                if(posteriors[j] > posteriors[maxPosteriorIndex]){	
                    maxPosteriorIndex = j;
                }
                
            }
            
            maxIndex[i] = getMinRiskIndexKmer(queryProfile, indices[i], posteriors);
            
            maxIndex[i] = maxPosteriorIndex;
            bestPosterior[i] = posteriors[maxIndex[i]];	
        }
        
        int saneDepth = sanityCheck(indices, maxIndex);
        
        simpleTax = "";
        int savedspot = 1;
        taxonProbabilityString = "";
        for(int i=1;i<=saneDepth;i++){
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            int confidenceScore = (int) (bestPosterior[i] * 100);
            if (confidenceScore >= confidenceThreshold) {
                if(indices[i][maxIndex[i]] != -1){
                    taxonProbabilityString += tree[indices[i][maxIndex[i]]]->getName() + "(" + toString(confidenceScore) + ");";
                    simpleTax += tree[indices[i][maxIndex[i]]]->getName() + ";";
                }
                else{
                    taxonProbabilityString += "unclassified(" + toString(confidenceScore) + ");";
                    simpleTax += "unclassified;";
                }
            }else { break; }
            savedspot = i;
        }
        
        
        
        for(int i=savedspot+1;i<numLevels;i++){
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            taxonProbabilityString += "unclassified(0);";
            simpleTax += "unclassified;";
        }
        
        return taxonProbabilityString;
	}
	catch(exception& e) {
		m->errorOut(e, "KmerTree", "getTaxonomy");
		exit(1);
	}
}


/**************************************************************************************************/

