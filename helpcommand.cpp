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
		cout << "The read.dist command parameter options are phylipfile or columnfile, namefile, cutoff and precision" << "\n";
		cout << "The read.dist command should be in the following format: " << "\n";
		cout << "read.dist(phylipfile=yourDistFile, namefile=yourNameFile, cutoff=yourCutoff, precision=yourPrecision) " << "\n";
		cout << "The phylipfile or columnfile parameter is required, but only one may be used.  If you use a columnfile the namefile is required. " << "\n";
		cout << "If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. distfile), '=' and parameters (i.e.yourDistfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "read.otu") {
		cout << "The read.otu command parameter options are listfile, rabundfile, sabundfile or orderfile." << "\n";
		cout << "The read.otu command should be in the following format: " << "\n";
		cout << "read.otu(listfile=yourListFile, orderfile=yourOrderFile) " << "\n";
		cout << "The read.otu requires one of hte following parameters: listfile, rabundfile or sabundfile. Only one may be used at a time." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "read.list") {
		cout << "The read.list command parameter options are listfile and groupfile." << "\n";
		cout << "The read.list command should be in the following format: " << "\n";
		cout << "read.list(listfile=yourListFile, groupfile=yourGroupFile) " << "\n";
		cout << "The listfile parameter and groupfile paramaters are required." << "\n";
		cout << "The read.list command parses a list file and separates it into groups." << "\n";
		cout << "It outputs a .shared file containing the otu information for each group as well as a .list file for each group." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "cluster") {
		cout << "The cluster command can only be executed after a successful read.phylip or read.column command." << "\n";
		cout << "The cluster command parameter options are method, cuttoff and precision. No parameters are required." << "\n";
		cout << "The cluster command should be in the following format: " << "\n";
		cout << "cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) " << "\n";
		cout << "The acceptable cluster methods are furthest, nearest and average.  If no method is provided then furthest is assumed." << "\n" << "\n";
	}else if (globaldata->helpRequest == "collect.single") {
		cout << "The collect.single command can only be executed after a successful read.list read.rabund or rad.sabund command. WITH ONE EXECEPTION. " << "\n";
		cout << "The collect.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The collect.single command parameters are label, line, freq, single.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.single command should be in the following format: " << "\n";
		cout << "collect.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, single=yourEstimators)." << "\n";
		cout << "Example collect(label=unique-.01-.03, line=0,5,10, iters=10000, freq=10, single=collect-chao-ace-jack)." << "\n";
		cout << "The default values for freq is 100, and single are sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "collect.shared") {
		cout << "The collect.shared command can only be executed after a successful read.shared command." << "\n";
		cout << "The collect.shared command parameters are label, line, freq, jumble, shared.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The collect.shared command should be in the following format: " << "\n";
		cout << "collect.shared(label=yourLabel, line=yourLines, freq=yourFreq, jumble=yourJumble, shared=yourEstimators)." << "\n";
		cout << "Example collect.shared(label=unique-.01-.03, line=0,5,10, freq=10, jumble=1, shared=sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN)." << "\n";
		cout << "The default values for jumble is 0 (meaning don’t jumble, if it’s set to 1 then it will jumble), freq is 100 and shared are sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "rarefaction.single") {
		cout << "The rarefaction.single command can only be executed after a successful read.list, read.rabund or read.sabund. WTIH ONE EXECEPTION." << "\n";
		cout << "The rarefaction.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The rarefaction.single command parameters are label, line, iters, freq, rarefaction.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction.single command should be in the following format: " << "\n";
		cout << "rarefaction.single(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, rarefaction=yourEstimators)." << "\n";
		cout << "Example rarefaction.single(label=unique-.01-.03, line=0,5,10, iters=10000, freq=10, rarefaction=rarefaction-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson)." << "\n";
		cout << "The default values for iters is 1000, freq is 100, and rarefaction is rarefaction which calculates the rarefaction curve for the observed richness." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "rarefaction.shared") {
		cout << "The rarefaction.shared command can only be executed after a successful read.shared command." << "\n";
		cout << "The rarefaction.shared command parameters are label, line, iters, jumble and sharedrarefaction.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The rarefaction command should be in the following format: " << "\n";
		cout << "rarefaction.shared(label=yourLabel, line=yourLines, iters=yourIters, jumble=yourJumble, sharedrarefaction=yourEstimators)." << "\n";
		cout << "Example rarefaction.shared(label=unique-.01-.03, line=0,5,10, iters=10000, jumble=1, sharedrarefaction=sharedobserved)." << "\n";
		cout << "The default values for iters is 1000, jumble is 0 (meaning don’t jumble, if it’s set to 1 then it will jumble), freq is 100, and sharedrarefaction is sharedobserved which calculates the shared rarefaction curve for the observed richness." << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.single") { 
		cout << "The summary.single command can only be executed after a successful read.list, read.rabund or read.sabund. WTIH ONE EXECEPTION." << "\n";
		cout << "The summary.single command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster." << "\n";
		cout << "The summary.single command parameters are label, line, summary.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.single command should be in the following format: " << "\n";
		cout << "summary.single(label=yourLabel, line=yourLines, summary=yourEstimators)." << "\n";
		cout << "Example summary.single(label=unique-.01-.03, line=0,5,10, summary=sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson)." << "\n";
		cout << "The default value summary is sobs-chao-ace-jack-bootstrap-shannon-npshannon-simpson" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "summary.shared") { 
		cout << "The summary.shared command can only be executed after a successful read.shared command." << "\n";
		cout << "The summary.shared command parameters are label, line, jumble and sharedsummary.  No parameters are required, but you may not use " << "\n";
		cout << "both the line and label parameters at the same time. The summary.shared command should be in the following format: " << "\n";
		cout << "summary.shared(label=yourLabel, line=yourLines, jumble=yourJumble, sharedsummary=yourEstimators)." << "\n";
		cout << "Example summary.shared(label=unique-.01-.03, line=0,5,10, jumble=1, sharedsummary=sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN)." << "\n";
		cout << "The default value for jumble is 0 (meaning don’t jumble, if it’s set to 1 then it will jumble) and sharedsummary is sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN" << "\n";
		cout << "The label and line parameters are used to analyze specific lines in your input." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. listfile), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";
	}else if (globaldata->helpRequest == "quit") {
		cout << "The quit command will terminate Dotur and should be in the following format: " << "\n";
		cout << "quit()" << "\n" << "\n";
	}else if (globaldata->helpRequest == "") {
		cout << "Valid commands are read.dist(), read.list(), read.otu(), cluster(), collect.single(), rarefaction.single(), summary.single(), collect.shared(), rarefaction.shared(), summary.shared(), quit(), help()." << "\n";
		cout << "For more information about a specific command type 'help(commandName)' i.e. 'help(read.phylip)'" << endl;
	}else {
		cout << "not a valid command" << endl;
	}
	
	cout << endl << "For further assistance please refer to the Mothur manual, or contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
	return 0;
}

//**********************************************************************************************************************/
