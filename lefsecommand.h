//
//  lefsecommand.h
//  Mothur
//
//  Created by SarahsWork on 6/12/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__lefsecommand__
#define __Mothur__lefsecommand__

#include "command.hpp"

/* 
 Columns = groups, rows are OTUs, class = design
 
 From http://huttenhower.sph.harvard.edu/galaxy/root?tool_id=lefse_upload
 Input data consist of a collection of m samples (columns) each made up of n numerical features (rows, typically normalized per-sample, red representing high values and green low). These samples are labeled with a class (taking two or more possible values) that represents the main biological hypothesis under investigation; they may also have one or more subclass labels reflecting within-class groupings.
 
 Step 1: the Kruskall-Wallis test analyzes all features, testing whether the values in different classes are differentially distributed. Features violating the null hypothesis are further analyzed in Step 2.
 Step 2: the pairwise Wilcoxon test checks whether all pairwise comparisons between subclasses within different classes significantly agree with the class level trend.
 Step 3: the resulting subset of vectors is used to build a Linear Discriminant Analysis model from which the relative difference among classes is used to rank the features. The final output thus consists of a list of features that are discriminative with respect to the classes, consistent with the subclass grouping within classes, and ranked according to the effect size with which they differentiate classes.
*/


#include "command.hpp"
#include "inputdata.h"
#include "designmap.h"

/**************************************************************************************************/

class LefseCommand : public Command {
public:
    LefseCommand(string);
    LefseCommand();
    ~LefseCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "lefse";			}
    string getCommandCategory()		{ return "OTU-Based Approaches";		}
    
    string getOutputPattern(string);
	string getHelpString();
    string getCitation() { return "http://www.mothur.org/wiki/Lefse"; }
    string getDescription()		{ return "brief description"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, allLines, wilc, wilcsamename, curv, subject, normMillion;
    string outputDir, sharedfile, designfile, mclass, subclass, classes, rankTec, multiClassStrat;
    vector<string> outputNames;
    set<string> labels;
    double anovaAlpha, wilcoxonAlpha, fBoots, ldaThreshold;
    int nlogs, iters, strict, minC;
    
    int process(vector<SharedRAbundFloatVector*>&, DesignMap&);
    int normalize(vector<SharedRAbundFloatVector*>&);
    map<int, double> runKruskalWallis(vector<SharedRAbundFloatVector*>&, DesignMap&);
    map<int, double> runWilcoxon(vector<SharedRAbundFloatVector*>&, DesignMap&, map<int, double>, map<string, set<string> >& class2SubClasses, map<string, vector<int> >& subClass2GroupIndex, map<string, string>);
    bool testOTUWilcoxon(map<string, set<string> >& class2SubClasses, vector<float> abunds, map<string, vector<int> >& subClass2GroupIndex, map<string, string>);
    map<int, double> testLDA(vector<SharedRAbundFloatVector*>&, map<int, double>, map<string, vector<int> >& class2GroupIndex, map<string, vector<int> >&);
    bool contastWithinClassesOrFewPerClass(vector< vector<double> >&, vector<int> rands, int minCl, map<string, vector<int> > class2GroupIndex,  map<int, string> indexToClass);
    vector< vector<double> > lda(vector< vector<double> >& adjustedLookup, vector<int> rand_s, map<int, string>& indexToClass, vector<string>);
    vector< vector<double> > getMeans(vector<SharedRAbundFloatVector*>& lookup, map<string, vector<int> >& class2GroupIndex);
    int printResults(vector< vector<double> >, map<int, double>, map<int, double>, string, vector<string>);
    
    //for testing
    bool printToCoutForRTesting(vector< vector<double> >& adjustedLookup, vector<int> rand_s, map<string, vector<int> >& class2GroupIndex, map<int, double> bins, map<string, vector<int> >&, vector<string>);
    int makeShared(int);
};

/**************************************************************************************************/




#endif /* defined(__Mothur__lefsecommand__) */
