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

HelpCommand::HelpCommand(){}

//**********************************************************************************************************************

HelpCommand::~HelpCommand(){}

//**********************************************************************************************************************

int HelpCommand::execute(){

	globaldata = GlobalData::getInstance();
	
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
		cout << "The collect.single command parameters are label, line, freq, single.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.single command should be in the following format: " << "\n";
		cout << "collect.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, single=yourEstimators)." << "\n";
		cout << "Example collect(label=unique-.01-.03, line=0,5,10, iters=10000, freq=10, single=collect-chao-ace-jack)." << "\n";
		cout << "The default values for freq is 100, and single are sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "collect.shared") {
		cout << "The collect.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The collect.shared command parameters are label, line, freq, jumble, shared.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.shared command should be in the following format: " << "\n";
		cout << "collect.shared(label=yourLabel, line=yourLines, freq=yourFreq, jumble=yourJumble, shared=yourEstimators)." << "\n";
		cout << "Example collect.shared(label=unique-.01-.03, line=0,5,10, freq=10, jumble=1, shared=sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN)." << "\n";
		cout << "The default values for jumble is 1 (meaning jumble, if it’s set to 0 then it will not jumble), freq is 100 and shared are sharedsobs-sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "rarefaction.single") {
		cout << "The rarefaction.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION." << "\n";
		cout << "The rarefaction.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The rarefaction.single command parameters are label, line, iters, freq, rarefaction.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction.single command should be in the following format: " << "\n";
		cout << "rarefaction.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, rarefaction=yourEstimators)." << "\n";
		cout << "Example rarefaction.single(label=unique-.01-.03, line=0,5,10, iters=10000, freq=10, rarefaction=sobs-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson)." << "\n";
		cout << "The default values for iters is 1000, freq is 100, and rarefaction is rarefaction which calculates the rarefaction curve for the observed richness." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "rarefaction.shared") {
		cout << "The rarefaction.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The rarefaction.shared command parameters are label, line, iters, jumble and sharedrarefaction.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction command should be in the following format: " << "\n";
		cout << "rarefaction.shared(label=yourLabel, line=yourLines, iters=yourIters, jumble=yourJumble, sharedrarefaction=yourEstimators)." << "\n";
		cout << "Example rarefaction.shared(label=unique-.01-.03, line=0,5,10, iters=10000, jumble=1, sharedrarefaction=sharedobserved)." << "\n";
		cout << "The default values for iters is 1000, jumble is 1 (meaning jumble, if it’s set to 0 then it will not jumble), freq is 100, and sharedrarefaction is sharedobserved which calculates the shared rarefaction curve for the observed richness." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.single") { 
		cout << "The summary.single command can only be executed after a successful read.otu WTIH ONE EXECEPTION." << "\n";
		cout << "The summary.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The summary.single command parameters are label, line, summary.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.single command should be in the following format: " << "\n";
		cout << "summary.single(label=yourLabel, line=yourLines, summary=yourEstimators)." << "\n";
		cout << "Example summary.single(label=unique-.01-.03, line=0,5,10, summary=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson)." << "\n";
		cout << "The default value summary is sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.shared") { 
		cout << "The summary.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The summary.shared command parameters are label, line, jumble and sharedsummary.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.shared command should be in the following format: " << "\n";
		cout << "summary.shared(label=yourLabel, line=yourLines, jumble=yourJumble, sharedsummary=yourEstimators)." << "\n";
		cout << "Example summary.shared(label=unique-.01-.03, line=0,5,10, jumble=1, sharedsummary=sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN)." << "\n";
		cout << "The default value for jumble is 1 (meaning jumble, if it’s set to 0 then it will not jumble) and sharedsummary is sharedsobs-sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "quit") {
		cout << "The quit command will terminate Dotur and should be in the following format: " << "\n";
		cout << "quit()" << "\n" << "\n";
	}else if (globaldata->helpRequest == "") {
		cout << "Valid commands are read.dist(), read.otu(), cluster(), deconvolute(), collect.single(), rarefaction.single(), summary.single(), collect.shared(), rarefaction.shared(), summary.shared(), quit(), help()." << "\n";
		cout << "For more information about a specific command type 'help(commandName)' i.e. 'help(read.dist)'" << endl;
	}else {
		cout << globaldata->helpRequest << " is not a valid command" << endl;
	}
	
	cout << endl << "For further assistance please refer to the Mothur manual, or contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
	return 0;
}

//**********************************************************************************************************************/
