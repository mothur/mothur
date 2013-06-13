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
 Columns = groups, rows are OTUs, class = design??
 
 From http://huttenhower.sph.harvard.edu/galaxy/root?tool_id=lefse_upload
 Input data consist of a collection of m samples (columns) each made up of n numerical features (rows, typically normalized per-sample, red representing high values and green low). These samples are labeled with a class (taking two or more possible values) that represents the main biological hypothesis under investigation; they may also have one or more subclass labels reflecting within-class groupings.
 
 Step 1: the Kruskall-Wallis test analyzes all features, testing whether the values in different classes are differentially distributed. Features violating the null hypothesis are further analyzed in Step 2.
 Step 2: the pairwise Wilcoxon test checks whether all pairwise comparisons between subclasses within different classes significantly agree with the class level trend.
 Step 3: the resulting subset of vectors is used to build a Linear Discriminant Analysis model from which the relative difference among classes is used to rank the features. The final output thus consists of a list of features that are discriminative with respect to the classes, consistent with the subclass grouping within classes, and ranked according to the effect size with which they differentiate classes.
*/

#endif /* defined(__Mothur__lefsecommand__) */
