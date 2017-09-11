//
//  alignTree.cpp
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#include "alignnode.h"
#include "aligntree.h"

/**************************************************************************************************/

AlignTree::AlignTree(string referenceFileName, string taxonomyFileName, int cutoff) : Classify(), confidenceThreshold(cutoff){
	try {
        AlignNode* newNode = new AlignNode("Root", 0);
        tree.push_back(newNode);			//	the tree is stored as a vector of elements of type TaxonomyNode
        
        string refTaxonomy;
        
        readTaxonomy(taxonomyFileName);
     
        ifstream referenceFile;
        m->openInputFile(referenceFileName, referenceFile);
        bool error = false;
        map<int, int> lengths;
        while(!referenceFile.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(referenceFile);  m->gobble(referenceFile);
            
            if (seq.getName() != "") {
                map<string, string>::iterator it = taxonomy.find(seq.getName());
                
                if (it != taxonomy.end()) {
                    refTaxonomy = it->second;		//	lookup the taxonomy string for the current reference sequence
                    string aligned = seq.getAligned();
                    lengths[aligned.length()] = 1;
                    if (lengths.size() > 1) { error = true; m->mothurOut("[ERROR]: reference sequences must be aligned to use the align method, quitting.\n"); break; }
                    addTaxonomyToTree(seq.getName(), refTaxonomy, aligned);
                }else {
                    m->mothurOut(seq.getName() + " is in your reference file, but not in your taxonomy file, please correct.\n"); error = true;
                }
            }
        }
        referenceFile.close();
        
        length = (lengths.begin())->first;  
           
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
            tree[i]->setTotalSeqs(dbSize);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "AlignTree", "AlignTree");
        exit(1);
    }
}

/**************************************************************************************************/

AlignTree::~AlignTree(){
	try {
        for(int i=0;i<tree.size();i++){
            delete tree[i];
        }
	}
    catch(exception& e) {
        m->errorOut(e, "AlignTree", "~AlignTree");
        exit(1);
    }
}	

/**************************************************************************************************/

int AlignTree::addTaxonomyToTree(string seqName, string& taxonomy, string& sequence){
	try {
        AlignNode* newNode;
        string taxonName = "";
        int treePosition = 0;							//	the root is element 0
        
        int level = 1;
        
        for(int i=0;i<taxonomy.length();i++){			//	step through taxonomy string...
            
            if (m->getControl_pressed()) { break; }
            
            if(taxonomy[i] == ';'){						//	looking for semicolons...
                
                if (taxonName == "") {  m->mothurOut(seqName + " has an error in the taxonomy.  This may be due to a ;;"); m->mothurOutEndLine(); m->setControl_pressed(true); }
                
                int newIndex = tree[treePosition]->getChildIndex(taxonName);	//	look to see if your current node already
                //	has a child with the new taxonName
                if(newIndex != -1)	{	treePosition = newIndex;	}		//	if you've seen it before, jump to that
                else {														//	 position in the tree
                    int newChildIndex = (int)tree.size();						//	otherwise, we'll have to create one...
                    tree[treePosition]->makeChild(taxonName, newChildIndex);
                    
                    newNode = new AlignNode(taxonName, level);
                    
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
        m->errorOut(e, "AlignTree", "addTaxonomyToTree");
        exit(1);
    }
}

/**************************************************************************************************/

int AlignTree::aggregateThetas(){
	try {
        vector<vector<int> > levelMatrix(numLevels+1);
        
        for(int i=0;i<tree.size();i++){
            if (m->getControl_pressed()) { return 0; }
            levelMatrix[tree[i]->getLevel()].push_back(i);
        }
		
        for(int i=numLevels-1;i>0;i--){
            if (m->getControl_pressed()) { return 0; }
            for(int j=0;j<levelMatrix[i].size();j++){
                
                AlignNode* holder = tree[levelMatrix[i][j]];
                
                tree[holder->getParent()]->addThetas(holder->getTheta(), holder->getNumSeqs());				
            }
        }
	    return 0;
	}
    catch(exception& e) {
        m->errorOut(e, "AlignTree", "aggregateThetas");
        exit(1);
    }
}

/**************************************************************************************************/

double AlignTree::getOutlierLogProbability(string& sequence){
	try {
        double count = 0;
        
        for(int i=0;i<sequence.length();i++){
            
            if(sequence[i] != '.'){	count++;	}
            
        }
        
        return count * log(0.2);
    }
    catch(exception& e) {
        m->errorOut(e, "AlignTree", "getOutlierLogProbability");
        exit(1);
    }
}

/**************************************************************************************************/

int AlignTree::getMinRiskIndexAlign(string& sequence, vector<int>& taxaIndices, vector<double>& probabilities){
	try {
        int numProbs = (int)probabilities.size();
        
        vector<double> G(numProbs, 0.2);	//a random sequence will, on average, be 20% similar to any other sequence
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
        m->errorOut(e, "AlignTree", "getMinRiskIndexAlign");
        exit(1);
    }

}

/**************************************************************************************************/

int AlignTree::sanityCheck(vector<vector<int> >& indices, vector<int>& maxIndices){
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
        m->errorOut(e, "AlignTree", "sanityCheck");
        exit(1);
    }
}

/**************************************************************************************************/

string AlignTree::getTaxonomy(Sequence* seq){
    try {
        string seqName = seq->getName(); string querySequence = seq->getAligned(); string taxonProbabilityString = "";
        if (querySequence.length() != length) {
            m->mothurOut("[ERROR]: " + seq->getName() + " has length " + toString(querySequence.length()) + ", reference sequences length is " + toString(length) + ". Are your sequences aligned? Sequences must be aligned to use the align search method.\n"); m->setControl_pressed(true); return "";
        }
        double logPOutlier = getOutlierLogProbability(querySequence);
        
        vector<vector<double> > pXgivenKj_D_j(numLevels);
        vector<vector<int> > indices(numLevels);
        for(int i=0;i<numLevels;i++){
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            pXgivenKj_D_j[i].push_back(logPOutlier);
            indices[i].push_back(-1);
        }
        
        
        for(int i=0;i<numTaxa;i++){
            //		cout << i << '\t' << tree[i]->getName() << '\t' << tree[i]->getLevel() << '\t' << tree[i]->getPxGivenkj_D_j(querySequence) << endl;
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            pXgivenKj_D_j[tree[i]->getLevel()].push_back(tree[i]->getPxGivenkj_D_j(querySequence));
            indices[tree[i]->getLevel()].push_back(i);
        }
        
        vector<double> sumLikelihood(numLevels, 0);
        vector<double> bestPosterior(numLevels, 0);
        vector<int> maxIndex(numLevels, 0);
        int maxPosteriorIndex;
        
        
        	//cout << "before best level" << endl;
        
        //let's find the best level and taxa within that level
        for(int i=0;i<numLevels;i++){ //go across all j's - from the root to genus
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            int numTaxaInLevel = (int)indices[i].size();
            
            		//cout << "numTaxaInLevel:\t" << numTaxaInLevel << endl;
            
            vector<double> posteriors(numTaxaInLevel, 0);		
            sumLikelihood[i] = getLogExpSum(pXgivenKj_D_j[i], maxPosteriorIndex);
            
            maxPosteriorIndex = 0;
            for(int j=0;j<numTaxaInLevel;j++){
                posteriors[j] = exp(pXgivenKj_D_j[i][j] - sumLikelihood[i]);
                
                if(posteriors[j] > posteriors[maxPosteriorIndex]){	
                    maxPosteriorIndex = j;
                }
                
            }
            
            maxIndex[i] = getMinRiskIndexAlign(querySequence, indices[i], posteriors);
            
            maxIndex[i] = maxPosteriorIndex;
            bestPosterior[i] = posteriors[maxIndex[i]];	
        }
        
        //	vector<double> pX_level(numLevels, 0);
        //	
        //	for(int i=0;i<numLevels;i++){
        //		pX_level[i] = pXgivenKj_D_j[i][maxIndex[i]] - tree[indices[i][maxIndex[i]]]->getNumSeqs();
        //	}
        //	
        //	int max_pLevel_X_index = -1;
        //	double pX_level_sum = getLogExpSum(pX_level, max_pLevel_X_index);
        //	double max_pLevel_X = exp(pX_level[max_pLevel_X_index] - pX_level_sum);
        //	
        //	vector<double> pLevel_X(numLevels, 0);
        //	for(int i=0;i<numLevels;i++){
        //		pLevel_X[i] = exp(pX_level[i] - pX_level_sum);
        //	}
        
        
        
        
        int saneDepth = sanityCheck(indices, maxIndex);
        
        simpleTax = "";
        int savedspot = 1;
        taxonProbabilityString = "";
        for(int i=1;i<=saneDepth;i++){
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            int confidenceScore = (int) (bestPosterior[i] * 100);
            if (confidenceScore >= confidenceThreshold) {
            if(indices[i][maxIndex[i]] != -1){
                taxonProbabilityString += tree[indices[i][maxIndex[i]]]->getName() + '(' + toString(confidenceScore) + ");";
                simpleTax += tree[indices[i][maxIndex[i]]]->getName() + ";";
                //			levelProbabilityOutput << tree[indices[i][maxIndex[i]]]->getName() << '(' << setprecision(6) << pLevel_X[i] << ");";
            }
            else{
                taxonProbabilityString + "unclassified" + '(' + toString(confidenceScore) + ");";
                //			levelProbabilityOutput << "unclassified" << '(' << setprecision(6) << pLevel_X[i] << ");";
                simpleTax += "unclassified;";
            }
            }else { break; }
            savedspot = i;
        }
        
        for(int i=savedspot+1;i<numLevels;i++){
            if (m->getControl_pressed()) { return taxonProbabilityString; }
            taxonProbabilityString + "unclassified(0);";
            simpleTax += "unclassified;";
        }
        
        return taxonProbabilityString;
	}
    catch(exception& e) {
        m->errorOut(e, "AlignTree", "getTaxonomy");
        exit(1);
    }
}


/**************************************************************************************************/
