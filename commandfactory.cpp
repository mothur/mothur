/*
 *  commandfactory.cpp
 *  
 *
 *  Created by Pat Schloss on 10/25/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "command.hpp"
#include "readdistcommand.h"
#include "readtreecommand.h"
#include "readotucommand.h"
#include "clustercommand.h"
#include "parselistcommand.h"
#include "collectcommand.h"
#include "collectsharedcommand.h"
#include "getgroupcommand.h"
#include "getlabelcommand.h"
#include "getlinecommand.h"
#include "rarefactcommand.h"
#include "summarycommand.h"
#include "summarysharedcommand.h"
#include "rarefactsharedcommand.h"
#include "quitcommand.h"
#include "helpcommand.h"
#include "commandfactory.hpp"
#include "deconvolutecommand.h"
#include "parsimonycommand.h"
#include "unifracunweightedcommand.h"
#include "unifracweightedcommand.h"
#include "libshuffcommand.h"
#include "heatmapcommand.h"
#include "heatmapsimcommand.h"
#include "filterseqscommand.h"
#include "mothur.h"
#include "venncommand.h"
#include "nocommands.h"
#include "binsequencecommand.h"
#include "getoturepcommand.h"
#include "treegroupscommand.h"
#include "bootstrapsharedcommand.h"
#include "concensuscommand.h"
#include "distancecommand.h"
#include "aligncommand.h"
#include "matrixoutputcommand.h"
#include "getsabundcommand.h"
#include "getrabundcommand.h"
#include "seqsummarycommand.h"
#include "screenseqscommand.h"
#include "reversecommand.h"
#include "trimseqscommand.h"

/***********************************************************/

/***********************************************************/
CommandFactory::CommandFactory(){
	string s = "";
	command = new NoCommand(s);
}
/***********************************************************/

/***********************************************************/
CommandFactory::~CommandFactory(){
	delete command;
}

/***********************************************************/

/***********************************************************/
//This function calls the appropriate command fucntions based on user input.
Command* CommandFactory::getCommand(string commandName, string optionString){
	try {
		delete command;   //delete the old command

		if(commandName == "read.dist")					{	command = new ReadDistCommand(optionString);			}
		else if(commandName == "read.otu")				{	command = new ReadOtuCommand(optionString);				}
		else if(commandName == "read.tree")				{	command = new ReadTreeCommand(optionString);			}
		else if(commandName == "cluster")				{	command = new ClusterCommand(optionString);				}
		else if(commandName == "unique.seqs")			{	command = new DeconvoluteCommand(optionString);			}
		else if(commandName == "parsimony")				{	command = new ParsimonyCommand(optionString);			}
		else if(commandName == "help")					{	command = new HelpCommand(optionString);				}
		else if(commandName == "quit")					{	command = new QuitCommand(optionString);				}
		else if(commandName == "collect.single")		{	command = new CollectCommand(optionString);				}
		else if(commandName == "collect.shared")		{	command = new CollectSharedCommand(optionString);		}
		else if(commandName == "rarefaction.single")	{	command = new RareFactCommand(optionString);			}
		else if(commandName == "rarefaction.shared")	{	command = new RareFactSharedCommand(optionString);		}
		else if(commandName == "summary.single")		{	command = new SummaryCommand(optionString);				}
		else if(commandName == "summary.shared")		{	command = new SummarySharedCommand(optionString);		}
		else if(commandName == "unifrac.weighted")		{	command = new UnifracWeightedCommand(optionString);		}
		else if(commandName == "unifrac.unweighted")	{	command = new UnifracUnweightedCommand(optionString);	}
		else if(commandName == "get.group")             {   command = new GetgroupCommand(optionString);			}
		else if(commandName == "get.label")             {   command = new GetlabelCommand(optionString);			}
		else if(commandName == "get.line")              {   command = new GetlineCommand(optionString);				}
		else if(commandName == "get.sabund")            {   command = new GetSAbundCommand(optionString);			}
		else if(commandName == "get.rabund")            {   command = new GetRAbundCommand(optionString);			}
		else if(commandName == "libshuff")              {   command = new LibShuffCommand(optionString);			}
		else if(commandName == "heatmap.bin")			{   command = new HeatMapCommand(optionString);				}
		else if(commandName == "heatmap.sim")			{   command = new HeatMapSimCommand(optionString);			}
		else if(commandName == "filter.seqs")			{   command = new FilterSeqsCommand(optionString);			}
		else if(commandName == "venn")					{   command = new VennCommand(optionString);				}
		else if(commandName == "bin.seqs")				{   command = new BinSeqCommand(optionString);				}
		else if(commandName == "get.oturep")			{   command = new GetOTURepCommand(optionString);			}
		else if(commandName == "tree.shared")			{   command = new TreeGroupCommand(optionString);			}
		else if(commandName == "dist.shared")			{   command = new MatrixOutputCommand(optionString);		}
		else if(commandName == "bootstrap.shared")		{   command = new BootSharedCommand(optionString);			}
		else if(commandName == "concensus")				{   command = new ConcensusCommand(optionString);			}
		else if(commandName == "dist.seqs")				{   command = new DistanceCommand(optionString);			}
		else if(commandName == "align.seqs")			{   command = new AlignCommand(optionString);				}
		else if(commandName == "summary.seqs")			{	command = new SeqSummaryCommand(optionString);			}
		else if(commandName == "screen.seqs")			{	command = new ScreenSeqsCommand(optionString);			}
		else if(commandName == "reverse.seqs")			{	command = new ReverseSeqsCommand(optionString);			}
		else if(commandName == "trim.seqs")				{	command = new TrimSeqsCommand(optionString);			}
		else											{	command = new NoCommand(optionString);					}

		return command;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CommandFactory class Function getCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CommandFactory class function getCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/***********************************************************/

