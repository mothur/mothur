/*
 *  helpcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "helpcommand.h"

//**********************************************************************************************************************

HelpCommand::HelpCommand(){
	globaldata = GlobalData::getInstance();
	validCommands = new ValidCommands();
	validCalcs = new ValidCalculators();
}

//**********************************************************************************************************************

HelpCommand::~HelpCommand(){}

//**********************************************************************************************************************

int HelpCommand::execute(){

	if (globaldata->helpRequest == "read.dist") {
		cout << "The read.dist command parameter options are phylip or column, group, name, cutoff and precision" << "\n";
		cout << "The read.dist command must be run before using the cluster or libshuff commands" << "\n";
		cout << "The read.dist command can be used in two ways.  The first is to read a phylip or column and run the cluster command" << "\n";
		cout << "For this use the read.dist command should be in the following format: " << "\n";
		cout << "read.dist(phylip=yourDistFile, name=yourNameFile, cutoff=yourCutoff, precision=yourPrecision) " << "\n";
		cout << "The phylip or column parameter is required, but only one may be used.  If you use a column file the name filename is required. " << "\n";
		cout << "If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed." << "\n";
		cout << "The second way to use the read.dist command is to read a phylip or column and a group, so you can use the libshuff command." << "\n";
		cout << "For this use the read.dist command should be in the following format: " << "\n";
		cout << "read.dist(phylip=yourPhylipfile, group=yourGroupFile). The cutoff and precision parameters are not valid with this use.  " << "\n";
		cout << "Note: No spaces between parameter labels (i.e. phylip), '=' and parameters (i.e.yourPhylipfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "read.otu") {
		cout << "The read.otu command must be run before you execute a collect.single, rarefaction.single, summary.single, " << "\n";
		cout << "collect.shared, rarefaction.shared or summary.shared command.   Mothur will generate a .list, .rabund and .sabund upon completion of the cluster command " << "\n";
		cout << "or you may use your own. The read.otu command parameter options are list, rabund, sabund, shared, group, order, line and label." << "\n";
		cout << "The read.otu command can be used in two ways.  The first is to read a list, rabund or sabund and run the collect.single, rarefaction.single or summary.single." << "\n";
		cout << "For this use the read.otu command should be in the following format: read.otu(list=yourListFile, order=yourOrderFile, label=yourLabels)." << "\n";
		cout << "The list, rabund or sabund parameter is required, but you may only use one of them." << "\n";
		cout << "The line and label parameters are optional but you may not use both the line and label parameters at the same time." << "\n";
		cout << "The label and line parameters are used to read specific lines in your input." << "\n";
		cout << "The second way to use the read.otu command is to read a list and a group, or a shared so you can use the collect.shared, rarefaction.shared or summary.shared commands." << "\n";
		cout << "In this case the read.otu command should be in the following format: read.otu(list=yourListFile, group=yourGroupFile, line=yourLines) or read.otu(shared=yourSharedFile).  " << "\n";
		cout << "The list parameter and group paramaters or the shared paremeter is required. When using the command the second way with a list and group file read.otu command parses the .list file" << "\n";
		cout << "and separates it into groups.  It outputs a .shared file containing the OTU information for each group. The read.otu command also outputs a .list file for each group. " << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "read.tree") {
		cout << "The read.tree command must be run before you execute a unifrac.weighted, unifrac.unweighted. " << "\n";
		cout << "It also must be run before using the parsimony command, unless you are using the randomtree parameter." << "\n";
		cout << "The read.tree command should be in the following format: read.tree(tree=yourTreeFile, group=yourGroupFile)." << "\n";
		cout << "The tree and group parameters are both required." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. tree), '=' and parameters (i.e.yourTreefile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "cluster") {
		cout << "The cluster command can only be executed after a successful read.dist command." << "\n";
		cout << "The cluster command parameter options are method, cuttoff and precision. No parameters are required." << "\n";
		cout << "The cluster command should be in the following format: " << "\n";
		cout << "cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) " << "\n";
		cout << "The acceptable cluster methods are furthest, nearest and average.  If no method is provided then furthest is assumed." << "\n" << "\n";
	}else if (globaldata->helpRequest == "deconvolute") {
		cout << "The deconvolute command reads a fastafile and creates a namesfile." << "\n";
		cout << "It creates a file where the first column is the groupname and the second column is a list of sequence names who have the same sequence. " << "\n";
		cout << "If the sequence is unique the second column will just contain its name. " << "\n";
		cout << "The deconvolute command parameter is fasta and it is required." << "\n";
		cout << "The deconvolute command should be in the following format: " << "\n";
		cout << "deconvolute(fasta=yourFastaFile) " << "\n";
	}else if (globaldata->helpRequest == "dist.seqs") {
		cout << "The dist.seqs command reads a file containing sequences and creates a distance file." << "\n";
		cout << "The dist.seqs command parameters are fasta, phylip, clustal, nexus, calc, countends, cutoff and processors.  " << "\n";
		cout << "You must use one of the following parameters for your filename: fasta, phylip, clustal or nexus. " << "\n";
		cout << "The calc parameter allows you to specify the method of calculating the distances.  Your options are: nogaps, onegap or eachgap. The default is onegap." << "\n";
		cout << "The countends parameter allows you to specify whether to include terminal gaps in distance.  Your options are: T or F. The default is T." << "\n";
		cout << "The cutoff parameter allows you to specify maximum distance to keep. The default is 1.0." << "\n";
		cout << "The processors parameter allows you to specify number of processors to use.  The default is 1, but you can use up to 4 processors." << "\n";
		cout << "The dist.seqs command should be in the following format: " << "\n";
		cout << "dist.seqs(fasta=yourFastaFile, calc=yourCalc, countends=yourEnds, cutoff= yourCutOff, processors=yourProcessors) " << "\n";
		cout << "Example dist.seqs(fasta=amazon.fasta, calc=eachgap, countends=F, cutoff= 2.0, processors=3)." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. calc), '=' and parameters (i.e.yourCalc)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "align.seqs") {
		cout << "The align.seqs command reads a file containing sequences and creates an alignment file and a report file." << "\n";
		cout << "The align.seqs command parameters are fasta, candidate, search, ksize, align, match, mismatch, gapopen and gapextend.  " << "\n";
		cout << "The template parameter is also required." << "\n";
		cout << "The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer." << "\n";
		cout << "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman." << "\n";
		cout << "The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 7." << "\n";
		cout << "The match parameter allows you to specify the bonus for having the same base. The default is 1.0." << "\n";
		cout << "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0." << "\n";
		cout << "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -1.0." << "\n";
		cout << "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -2.0." << "\n";
		cout << "The align.seqs command should be in the following format: " << "\n";
		cout << "align.seqs(fasta=yourTemplateFile, candidate=yourCandidateFile, align=yourAlignmentMethod, search=yourSearchmethod, ksize=yourKmerSize, match=yourMatchBonus, mismatch=yourMismatchpenalty, gapopen=yourGapopenPenalty, gapextend=yourGapExtendPenalty) " << "\n";
		cout << "Example align.seqs(candidate=candidate.fasta, fasta=core.filtered, align=kmer, search=gotoh, ksize=8, match=2.0, mismatch=3.0, gapopen=-2.0, gapextend=-1.0)" << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "collect.single") {
		cout << "The collect.single command can only be executed after a successful read.otu command. WITH ONE EXECEPTION. " << "\n";
		cout << "The collect.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The collect.single command parameters are label, line, freq, calc and abund.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.single command should be in the following format: " << "\n";
		cout << "collect.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, calc=yourEstimators)." << "\n";
		cout << "Example collect(label=unique-.01-.03, line=0-5-10, iters=10000, freq=10, calc=sobs-chao-ace-jack)." << "\n";
		cout << "The default values for freq is 100, and calc are sobs-chao-ace-jack-shannon-npshannon-simpson." << "\n";
		validCalcs->printCalc("single", cout);
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "collect.shared") {
		cout << "The collect.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The collect.shared command parameters are label, line, freq, calc and groups.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.shared command should be in the following format: " << "\n";
		cout << "collect.shared(label=yourLabel, line=yourLines, freq=yourFreq, calc=yourEstimators, groups=yourGroups)." << "\n";
		cout << "Example collect.shared(label=unique-.01-.03, line=0-5-10, freq=10, groups=B-C, calc=sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan)." << "\n";
		cout << "The default values for freq is 100 and calc are sharedsobs-sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan." << "\n";
		cout << "The default value for groups is all the groups in your groupfile." << "\n";
		validCalcs->printCalc("shared", cout);
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "get.group") {
		cout << "The get.group command can only be executed after a successful read.otu command of a group file." << "\n";
		cout << "You may not use any parameters with the get.group command." << "\n";
		cout << "The get.group command should be in the following format: " << "\n";
		cout << "get.group()" << "\n";
		cout << "Example get.group()." << "\n";
	}else if (globaldata->helpRequest == "get.label") {
		cout << "The get.label command can only be executed after a successful read.otu command." << "\n";
		cout << "You may not use any parameters with the get.label command." << "\n";
		cout << "The get.label command should be in the following format: " << "\n";
		cout << "get.label()" << "\n";
		cout << "Example get.label()." << "\n";
	}else if (globaldata->helpRequest == "get.line") {
		cout << "The get.line command can only be executed after a successful read.otu command." << "\n";
		cout << "You may not use any parameters with the get.line command." << "\n";
		cout << "The get.line command should be in the following format: " << "\n";
		cout << "get.line()" << "\n";
		cout << "Example get.line()." << "\n";
	}else if (globaldata->helpRequest == "rarefaction.single") {
		cout << "The rarefaction.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION." << "\n";
		cout << "The rarefaction.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The rarefaction.single command parameters are label, line, iters, freq, calc and abund.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction.single command should be in the following format: " << "\n";
		cout << "rarefaction.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, calc=yourEstimators)." << "\n";
		cout << "Example rarefaction.single(label=unique-.01-.03, line=0-5-10, iters=10000, freq=10, calc=sobs-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson)." << "\n";
		cout << "The default values for iters is 1000, freq is 100, and calc is rarefaction which calculates the rarefaction curve for the observed richness." << "\n";
		validCalcs->printCalc("rarefaction", cout);
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "rarefaction.shared") {
		cout << "The rarefaction.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The rarefaction.shared command parameters are label, line, iters, jumble, groups and calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction command should be in the following format: " << "\n";
		cout << "rarefaction.shared(label=yourLabel, line=yourLines, iters=yourIters, jumble=yourJumble, calc=yourEstimators, groups=yourGroups)." << "\n";
		cout << "Example rarefaction.shared(label=unique-.01-.03, line=0-5-10, iters=10000, jumble=1, groups=B-C, calc=sharedobserved)." << "\n";
		cout << "The default values for iters is 1000, jumble is 1 (meaning jumble, if it’s set to 0 then it will not jumble), freq is 100, and calc is sharedobserved which calculates the shared rarefaction curve for the observed richness." << "\n";
		cout << "The default value for groups is all the groups in your groupfile." << "\n";
		validCalcs->printCalc("sharedrarefaction", cout);
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.single") { 
		cout << "The summary.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION." << "\n";
		cout << "The summary.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The summary.single command parameters are label, line, calc, abund.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.single command should be in the following format: " << "\n";
		cout << "summary.single(label=yourLabel, line=yourLines, calc=yourEstimators)." << "\n";
		cout << "Example summary.single(label=unique-.01-.03, line=0,5,10, calc=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson)." << "\n";
		validCalcs->printCalc("summary", cout);
		cout << "The default value calc is sobs-chao-ace-jack-shannon-npshannon-simpson" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. line), '=' and parameters (i.e.yourLines)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.shared") { 
		cout << "The summary.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The summary.shared command parameters are label, line and calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.shared command should be in the following format: " << "\n";
		cout << "summary.shared(label=yourLabel, line=yourLines, calc=yourEstimators, groups=yourGroups)." << "\n";
		cout << "Example summary.shared(label=unique-.01-.03, line=0,5,10, groups=B-C, calc=sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan)." << "\n";
		validCalcs->printCalc("sharedsummary", cout);
		cout << "The default value for calc is sharedsobs-sharedchao-sharedace-jabund-sorensonabund-jclass-sorclass-jest-sorest-thetayc-thetan" << "\n";
		cout << "The default value for groups is all the groups in your groupfile." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. line), '=' and parameters (i.e.yourLines)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "parsimony") { 
		cout << "The parsimony command can only be executed after a successful read.tree command, unless you use the random parameter." << "\n";
		cout << "The parsimony command parameters are random, groups and iters.  No parameters are required." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree." << "\n";
		cout << "The parsimony command should be in the following format: parsimony(random=yourOutputFilename, groups=yourGroups, iters=yourIters)." << "\n";
		cout << "Example parsimony(random=out, iters=500)." << "\n";
		cout << "The default value for random is "" (meaning you want to use the trees in your inputfile, randomtree=out means you just want the random distribution of trees outputted to out.rd_parsimony)," << "\n";
		cout << "and iters is 1000.  The parsimony command output two files: .parsimony and .psummary their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. random), '=' and parameters (i.e.yourOutputFilename)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "unifrac.weighted") { 
		cout << "The unifrac.weighted command can only be executed after a successful read.tree command." << "\n";
		cout << "The unifrac.weighted command parameters are groups and iters.  No parameters are required." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree." << "\n";
		cout << "The unifrac.weighted command should be in the following format: unifrac.weighted(groups=yourGroups, iters=yourIters)." << "\n";
		cout << "Example unifrac.weighted(groups=A-B-C, iters=500)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, and iters is 1000." << "\n";
		cout << "The unifrac.weighted command output two files: .weighted and .wsummary their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "unifrac.unweighted") { 
		cout << "The unifrac.unweighted command can only be executed after a successful read.tree command." << "\n";
		cout << "The unifrac.unweighted command parameters are groups and iters.  No parameters are required." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree." << "\n";
		cout << "The unifrac.unweighted command should be in the following format: unifrac.unweighted(groups=yourGroups, iters=yourIters)." << "\n";
		cout << "Example unifrac.unweighted(groups=A-B-C, iters=500)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, and iters is 1000." << "\n";
		cout << "The unifrac.unweighted command output two files: .unweighted and .uwsummary their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "libshuff") { 
		cout << "The libshuff command can only be executed after a successful read.dist command." << "\n";
		cout << "The libshuff command parameters are groups, iters, step, form and cutoff.  No parameters are required." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random matrices you would like compared to your matrix." << "\n";
		cout << "The step parameter allows you to specify change in distance you would like between each output if you are using the discrete form." << "\n";
		cout << "The form parameter allows you to specify if you would like to analyze your matrix using the discrete or integral form. Your options are integral or discrete." << "\n";
		cout << "The libshuff command should be in the following format: libshuff(groups=yourGroups, iters=yourIters, cutOff=yourCutOff, form=yourForm, step=yourStep)." << "\n";
		cout << "Example libshuff(groups=A-B-C, iters=500, form=discrete, step=0.01, cutOff=2.0)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, iters is 10000, cutoff is 1.0, form is integral and step is 0.01." << "\n";
		cout << "The libshuff command output two files: .coverage and .slsummary their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e.yourIters)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "heatmap") { 
		cout << "The heatmap command can only be executed after a successful read.otu command." << "\n";
		cout << "The heatmap command parameters are groups, sorted, scale, line and label.  No parameters are required, but you may not use line and label at the same time." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap." << "\n";
		cout << "The sorted parameter allows you to choose to see the file with the shared otus at the top or the otus in the order they appear in your input file. " << "\n";
		cout << "The scale parameter allows you to choose the range of color your bin information will be displayed with." << "\n";
		cout << "The group names are separated by dashes. The line and label allow you to select what distance levels you would like a heatmap created for, and are also separated by dashes." << "\n";
		cout << "The heatmap command should be in the following format: heatmap(groups=yourGroups, sorted=yourSorted, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example heatmap(groups=A-B-C, line=1-3-5, sorted=F, scale=log10)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, and all lines in your inputfile will be used." << "\n";
		cout << "The default value for sorted is T meaning you want the shared otus on top, you may change it to F meaning the exact representation of your input file." << "\n";
		cout << "The default value for scale is log10; your other options are log2 and linear." << "\n";
		cout << "The heatmap command outputs a .svg file for each line or label you specify." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "venn") { 
		cout << "The venn command can only be executed after a successful read.otu command." << "\n";
		cout << "The venn command parameters are groups, calc, line and label.  No parameters are required, but you may not use line and label at the same time." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like included in your venn diagram, you may only use a maximum of 4 groups." << "\n";
		cout << "The group names are separated by dashes. The line and label allow you to select what distance levels you would like a venn diagram created for, and are also separated by dashes." << "\n";
		cout << "The venn command should be in the following format: venn(groups=yourGroups, calc=yourCalcs, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example venn(groups=A-B-C, line=1-3-5, calc=sharedsobs-sharedchao)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile up to 4, and all lines in your inputfile will be used." << "\n";
		cout << "The default value for calc is sobs if you have only read a list file or if you have selected only one group, and sharedsobs if you have multiple groups." << "\n";
		cout << "The default available estimators for calc are sobs, chao and ace if you have only read a list file, and sharedsobs, sharedchao and sharedace if you have read a list and group file or a shared file." << "\n";
		cout << "The only estmiator available four 4 groups is sharedsobs." << "\n";
		cout << "The venn command outputs a .svg file for each calculator you specify at each distance you choose." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "tree.shared") { 
		cout << "The tree.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The tree.shared command parameters are groups, calc, line and label.  The calc parameter is required, and you may not use line and label at the same time." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like included used." << "\n";
		cout << "The group names are separated by dashes. The line and label allow you to select what distance levels you would like trees created for, and are also separated by dashes." << "\n";
		cout << "The tree.shared command should be in the following format: tree.shared(groups=yourGroups, calc=yourCalcs, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example tree.shared(groups=A-B-C, line=1-3-5, calc=jabund-sorabund)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile." << "\n";
		cout << "There is no default value for calc." << "\n";
		validCalcs->printCalc("treegroup", cout);
		cout << "The tree.shared command outputs a .tre file for each calculator you specify at each distance you choose." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "dist.shared") { 
		cout << "The dist.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The dist.shared command parameters are groups, calc, line and label.  The calc parameter is required, and you may not use line and label at the same time." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like included used." << "\n";
		cout << "The group names are separated by dashes. The line and label allow you to select what distance levels you would like distance matrices created for, and are also separated by dashes." << "\n";
		cout << "The dist.shared command should be in the following format: dist.shared(groups=yourGroups, calc=yourCalcs, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example dist.shared(groups=A-B-C, line=1-3-5, calc=jabund-sorabund)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile." << "\n";
		cout << "The default value for calc is jclass and thetayc." << "\n";
		validCalcs->printCalc("matrix", cout);
		cout << "The dist.shared command outputs a .dist file for each calculator you specify at each distance you choose." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "bootstrap.shared") { 
		cout << "The bootstrap.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The bootstrap.shared command parameters are groups, calc, iters, line and label.  The calc parameter is required, and you may not use line and label at the same time." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like included used." << "\n";
		cout << "The group names are separated by dashes. The line and label allow you to select what distance levels you would like trees created for, and are also separated by dashes." << "\n";
		cout << "The bootstrap.shared command should be in the following format: bootstrap.shared(groups=yourGroups, calc=yourCalcs, line=yourLines, label=yourLabels, iters=yourIters)." << "\n";
		cout << "Example bootstrap.shared(groups=A-B-C, line=1-3-5, calc=jabund-sorabund, iters=100)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile." << "\n";
		cout << "There is no default value for calc. The default for iters is 1000." << "\n";
		validCalcs->printCalc("boot", cout);
		cout << "The bootstrap.shared command outputs a .tre file for each calculator you specify at each distance you choose containing iters number of trees." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "concensus") { 
		cout << "The concensus command can only be executed after a successful read.tree command." << "\n";
		cout << "The concensus command has no parameters." << "\n";
		cout << "The concensus command should be in the following format: concensus()." << "\n";
		cout << "The concensus command output two files: .concensus.tre and .concensuspairs." << "\n";
		cout << "The .concensus.tre file contains the concensus tree of the trees in your input file." << "\n";
		cout << "The branch lengths are the percentage of trees in your input file that had the given pair." << "\n";
		cout << "The .concensuspairs file contains a list of the internal nodes in your tree.  For each node, the pair that was used in the concensus tree " << "\n";
		cout << "is reported with its percentage, as well as the other pairs that were seen for that node but not used and their percentages." << "\n" << "\n";
	}else if (globaldata->helpRequest == "bin.seqs") { 
		cout << "The bin.seqs command can only be executed after a successful read.otu command of a list file." << "\n";
		cout << "The bin.seqs command parameters are fasta, name, line, label and group.  The fasta parameter is required, and you may not use line and label at the same time." << "\n";
		cout << "The line and label allow you to select what distance levels you would like a output files created for, and are separated by dashes." << "\n";
		cout << "The bin.seqs command should be in the following format: bin.seqs(fasta=yourFastaFile, name=yourNamesFile, group=yourGroupFile, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example bin.seqs(fasta=amazon.fasta, group=amazon.groups, line=1-3-5, name=amazon.names)." << "\n";
		cout << "The default value for line and label are all lines in your inputfile." << "\n";
		cout << "The bin.seqs command outputs a .fasta file for each distance you specify appending the OTU number to each name." << "\n";
		cout << "If you provide a groupfile, then it also appends the sequences group to the name." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "get.oturep") { 
		cout << "The get.oturep command can only be executed after a successful read.dist command." << "\n";
		cout << "The get.oturep command parameters are list, fasta, name, group, line and label.  The fasta and list parameters are required, and you may not use line and label at the same time." << "\n";
		cout << "The line and label allow you to select what distance levels you would like a output files created for, and are separated by dashes." << "\n";
		cout << "The get.oturep command should be in the following format: get.oturep(fasta=yourFastaFile, list=yourListFile, name=yourNamesFile, group=yourGroupFile, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example get.oturep(fasta=amazon.fasta, list=amazon.fn.list, group=amazon.groups, line=1-3-5, name=amazon.names)." << "\n";
		cout << "The default value for line and label are all lines in your inputfile." << "\n";
		cout << "The get.oturep command outputs a .fastarep file for each distance you specify, selecting one OTU representative for each bin." << "\n";
		cout << "If you provide a groupfile, then it also appends the names of the groups present in that bin." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "quit") {
		cout << "The quit command will terminate mothur and should be in the following format: " << "\n";
		cout << "quit()" << "\n" << "\n";
	}else if (globaldata->helpRequest == "") {
		validCommands->printCommands(cout);
		cout << "For more information about a specific command type 'help(commandName)' i.e. 'help(read.dist)'" << endl;
	}else {
		cout << globaldata->helpRequest << " is not a valid command" << endl;
	}
	
	cout << endl << "For further assistance please refer to the Mothur manual on our wiki at http://schloss.micro.umass.edu/mothur/, or contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
	return 0;
}

//**********************************************************************************************************************/
