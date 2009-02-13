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
}

//**********************************************************************************************************************

HelpCommand::~HelpCommand(){}

//**********************************************************************************************************************

int HelpCommand::execute(){

	if (globaldata->helpRequest == "read.dist") {
		cout << "The read.dist command parameter options are phylip or column, name, cutoff and precision" << "\n";
		cout << "The read.dist command should be in the following format: " << "\n";
		cout << "read.dist(phylip=yourDistFile, name=yourNameFile, cutoff=yourCutoff, precision=yourPrecision) " << "\n";
		cout << "The phylip or column parameter is required, but only one may be used.  If you use a column file the name filename is required. " << "\n";
		cout << "If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. dist), '=' and parameters (i.e.yourDistfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "read.otu") {
		cout << "The read.otu command must be run before you execute a collect.single, rarefaction.single, summary.single, " << "\n";
		cout << "collect.shared, rarefaction.shared or summary.shared command.   Mothur will generate a .list, .rabund and .sabund upon completion of the cluster command " << "\n";
		cout << "or you may use your own. The read.otu command parameter options are list, rabund, sabund, group and order." << "\n";
		cout << "The read.otu command can be used in two ways.  The first is to read a list, rabund or sabund and run the collect.single, rarefaction.single or summary.single." << "\n";
		cout << "For this use the read.otu command should be in the following format: read.otu(list=yourListFile, order=yourOrderFile)." << "\n";
		cout << "The list, rabund or sabund parameter is required, but you may only use one of them." << "\n";
		cout << "The second way to use the read.otu command is to read a list and a group so you can use the collect.shared, rarefaction.shared or summary.shared commands." << "\n";
		cout << "In this case the read.otu command should be in the following format: read.otu(list=yourListFile, group=yourGroupFile).  " << "\n";
		cout << "The list parameter and group paramaters are required. When using the command the second way read.otu command parses the .list file" << "\n";
		cout << "and separates it into groups.  It outputs a .shared file containing the OTU information for each group. The read.otu command also outputs a .list file for each group. " << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "read.tree") {
		cout << "The read.tree command must be run before you execute a unifrac.weighted, unifrac.unweighted. " << "\n";
		cout << "It also must be run before using the parsimony command, unless you are using the randomtree parameter." << "\n";
		cout << "The read.tree command should be in the following format: read.tree(tree=yourTreeFile, group=yourGroupFile)." << "\n";
		cout << "The tree and group parameters are both required." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. tree), '=' and parameters (i.e.yourTreefile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "read.shared") {
		cout << "The read.shared must be run before you execute a collect.shared, rarefaction.shared or summary.shared command." << "\n";
		cout << "The read.shared command is used to read a shared file. The read.shared should be entered in the following format:" << "\n";
		cout << "read.shared(shared=yourSharedFile). The shared parameter is required." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. shared), '=' and parameters (i.e.yourSharedfile)." << "\n" << "\n";
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
	}else if (globaldata->helpRequest == "collect.single") {
		cout << "The collect.single command can only be executed after a successful read.otu command. WITH ONE EXECEPTION. " << "\n";
		cout << "The collect.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The collect.single command parameters are label, line, freq, calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.single command should be in the following format: " << "\n";
		cout << "collect.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, calc=yourEstimators)." << "\n";
		cout << "Example collect(label=unique-.01-.03, line=0,5,10, iters=10000, freq=10, calc=sobs-chao-ace-jack)." << "\n";
		cout << "The default values for freq is 100, and calc are sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "collect.shared") {
		cout << "The collect.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The collect.shared command parameters are label, line, freq, jumble, calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.shared command should be in the following format: " << "\n";
		cout << "collect.shared(label=yourLabel, line=yourLines, freq=yourFreq, jumble=yourJumble, calc=yourEstimators)." << "\n";
		cout << "Example collect.shared(label=unique-.01-.03, line=0,5,10, freq=10, jumble=1, calc=sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN)." << "\n";
		cout << "The default values for jumble is 1 (meaning jumble, if it’s set to 0 then it will not jumble), freq is 100 and calc are sharedsobs-sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "rarefaction.single") {
		cout << "The rarefaction.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION." << "\n";
		cout << "The rarefaction.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The rarefaction.single command parameters are label, line, iters, freq, calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction.single command should be in the following format: " << "\n";
		cout << "rarefaction.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, calc=yourEstimators)." << "\n";
		cout << "Example rarefaction.single(label=unique-.01-.03, line=0,5,10, iters=10000, freq=10, calc=sobs-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson)." << "\n";
		cout << "The default values for iters is 1000, freq is 100, and calc is rarefaction which calculates the rarefaction curve for the observed richness." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "rarefaction.shared") {
		cout << "The rarefaction.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The rarefaction.shared command parameters are label, line, iters, jumble and calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction command should be in the following format: " << "\n";
		cout << "rarefaction.shared(label=yourLabel, line=yourLines, iters=yourIters, jumble=yourJumble, calc=yourEstimators)." << "\n";
		cout << "Example rarefaction.shared(label=unique-.01-.03, line=0,5,10, iters=10000, jumble=1, calc=sharedobserved)." << "\n";
		cout << "The default values for iters is 1000, jumble is 1 (meaning jumble, if it’s set to 0 then it will not jumble), freq is 100, and calc is sharedobserved which calculates the shared rarefaction curve for the observed richness." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.single") { 
		cout << "The summary.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION." << "\n";
		cout << "The summary.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The summary.single command parameters are label, line, calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.single command should be in the following format: " << "\n";
		cout << "summary.single(label=yourLabel, line=yourLines, calc=yourEstimators)." << "\n";
		cout << "Example summary.single(label=unique-.01-.03, line=0,5,10, calc=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson)." << "\n";
		cout << "The default value calc is sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. line), '=' and parameters (i.e.yourLines)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.shared") { 
		cout << "The summary.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The summary.shared command parameters are label, line, jumble and calc.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.shared command should be in the following format: " << "\n";
		cout << "summary.shared(label=yourLabel, line=yourLines, jumble=yourJumble, calc=yourEstimators)." << "\n";
		cout << "Example summary.shared(label=unique-.01-.03, line=0,5,10, jumble=1, calc=sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN)." << "\n";
		cout << "The default value for jumble is 1 (meaning jumble, if it’s set to 0 then it will not jumble) and calc is sharedsobs-sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. line), '=' and parameters (i.e.yourLines)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "parsimony") { 
		cout << "The parsimony command can only be executed after a successful read.tree command, unless you use the randomtree parameter." << "\n";
		cout << "The parsimony command parameters are randomtree and iters.  No parameters are required." << "\n";
		cout << "The parsimony command should be in the following format: parsimony(randomtree=yourRandomTreeValue, iters=yourIters)." << "\n";
		cout << "Example parsimony(randomtree=1, iters=500)." << "\n";
		cout << "The default value for randomTree is 0 (meaning you want to use the trees in your inputfile, randomtree=1 means you just want the random distribution of trees)," << "\n";
		cout << "and iters is 1000.  The parsimony command output three files: .parsimony, .psummary and .pdistrib, their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "unifrac.weighted") { 
		cout << "The unifrac.weighted command can only be executed after a successful read.tree command." << "\n";
		cout << "The unifrac.weighted command parameters are groups and iters.  No parameters are required." << "\n";
		cout << "The groups paramter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree." << "\n";
		cout << "The unifrac.weighted command should be in the following format: unifrac.weighted(groups=yourGroups, iters=yourIters)." << "\n";
		cout << "Example unifrac.weighted(groups=A-B-C, iters=500)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, and iters is 1000." << "\n";
		cout << "The unifrac.weighted command output three files: .weighted, .wsummary and .wdistrib, their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "unifrac.unweighted") { 
		cout << "The unifrac.unweighted command can only be executed after a successful read.tree command." << "\n";
		cout << "The unifrac.unweighted command parameters are groups and iters.  No parameters are required." << "\n";
		cout << "The groups paramter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree." << "\n";
		cout << "The unifrac.unweighted command should be in the following format: unifrac.unweighted(groups=yourGroups, iters=yourIters)." << "\n";
		cout << "Example unifrac.unweighted(groups=A-B-C, iters=500)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, and iters is 1000." << "\n";
		cout << "The unifrac.unweighted command output three files: .unweighted, .uwsummary and .uwdistrib, their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "quit") {
		cout << "The quit command will terminate Dotur and should be in the following format: " << "\n";
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
