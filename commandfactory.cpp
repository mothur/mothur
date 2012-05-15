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
#include "collectcommand.h"
#include "collectsharedcommand.h"
#include "getgroupcommand.h"
#include "getlabelcommand.h"
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
#include "venncommand.h"
#include "nocommands.h"
#include "binsequencecommand.h"
#include "getoturepcommand.h"
#include "treegroupscommand.h"
#include "distancecommand.h"
#include "aligncommand.h"
#include "matrixoutputcommand.h"
#include "getsabundcommand.h"
#include "getrabundcommand.h"
#include "seqsummarycommand.h"
#include "screenseqscommand.h"
#include "reversecommand.h"
#include "trimseqscommand.h"
#include "mergefilecommand.h"
#include "listseqscommand.h"
#include "getseqscommand.h"
#include "removeseqscommand.h"
#include "systemcommand.h"
#include "secondarystructurecommand.h"
#include "getsharedotucommand.h"
#include "getlistcountcommand.h"
#include "hclustercommand.h"
#include "classifyseqscommand.h"
#include "phylotypecommand.h"
#include "mgclustercommand.h"
#include "preclustercommand.h"
#include "pcoacommand.h"
#include "otuhierarchycommand.h"
#include "setdircommand.h"
#include "parselistscommand.h"
#include "chimeraccodecommand.h"
#include "chimeracheckcommand.h"
#include "chimeraslayercommand.h"
#include "chimerapintailcommand.h"
#include "chimerabellerophoncommand.h"
#include "chimerauchimecommand.h"
#include "setlogfilecommand.h"
#include "phylodiversitycommand.h"
#include "makegroupcommand.h"
#include "chopseqscommand.h"
#include "clearcutcommand.h"
#include "catchallcommand.h"
#include "splitabundcommand.h"
#include "clustersplitcommand.h"
#include "classifyotucommand.h"
#include "degapseqscommand.h"
#include "getrelabundcommand.h"
#include "sensspeccommand.h"
#include "sffinfocommand.h"
#include "seqerrorcommand.h"
#include "normalizesharedcommand.h"
#include "metastatscommand.h"
#include "splitgroupscommand.h"
#include "clusterfragmentscommand.h"
#include "getlineagecommand.h"
#include "removelineagecommand.h"
#include "parsefastaqcommand.h"
#include "pipelinepdscommand.h"
#include "deuniqueseqscommand.h"
#include "pairwiseseqscommand.h"
#include "clusterdoturcommand.h"
#include "subsamplecommand.h"
#include "removegroupscommand.h"
#include "getgroupscommand.h"
#include "getotuscommand.h"
#include "removeotuscommand.h"
#include "indicatorcommand.h"
#include "consensusseqscommand.h"
#include "trimflowscommand.h"
#include "corraxescommand.h"
#include "shhhercommand.h"
#include "pcacommand.h"
#include "nmdscommand.h"
#include "removerarecommand.h"
#include "mergegroupscommand.h"
#include "amovacommand.h"
#include "homovacommand.h"
#include "mantelcommand.h"
#include "makefastqcommand.h"
#include "anosimcommand.h"
#include "getcurrentcommand.h"
#include "setcurrentcommand.h"
#include "sharedcommand.h"
#include "getcommandinfocommand.h"
#include "deuniquetreecommand.h"
#include "countseqscommand.h"
#include "countgroupscommand.h"
#include "clearmemorycommand.h"
#include "summarytaxcommand.h"
#include "chimeraperseuscommand.h"
#include "shhhseqscommand.h"
#include "summaryqualcommand.h"
#include "otuassociationcommand.h"
#include "sortseqscommand.h"
#include "classifytreecommand.h"
#include "cooccurrencecommand.h"
#include "pcrseqscommand.h"
#include "createdatabasecommand.h"
#include "makebiomcommand.h"
#include "getcoremicrobiomecommand.h"
#include "listotulabelscommand.h"

/*******************************************************/

/******************************************************/
CommandFactory* CommandFactory::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new CommandFactory();
	}
	return _uniqueInstance;
}
/***********************************************************/

/***********************************************************/
//note: This class is resposible for knowing which commands are mpiEnabled,
//If a command is not enabled only process 0 will execute the command. 
//This avoids redundant outputs on pieces of code we have not paralellized. 
//If you add mpi code to a existing command you need to modify the list below or the code will hang on MPI blocking commands like FIle_open. 
//example:  commands["dist.seqs"] = "MPIEnabled";

CommandFactory::CommandFactory(){
	string s = "";
	m = MothurOut::getInstance();
	
	command = new NoCommand(s);
	shellcommand = new NoCommand(s);
	pipecommand = new NoCommand(s);
	
	outputDir = ""; inputDir = "";
	logFileName = "";
	append = false;
	
	//initialize list of valid commands
	commands["read.dist"]			= "read.dist"; 
	commands["read.otu"]			= "read.otu";
	commands["read.tree"]			= "read.tree"; 
	commands["make.shared"]			= "make.shared"; 
	commands["bin.seqs"]			= "bin.seqs"; 
	commands["get.oturep"]			= "get.oturep";
	commands["cluster"]				= "cluster"; 
	commands["unique.seqs"]			= "unique.seqs"; 
	commands["dist.shared"]			= "dist.shared";
	commands["collect.single"]		= "collect.single"; 
	commands["collect.shared"]		= "collect.shared"; 
	commands["rarefaction.single"]	= "rarefaction.single"; 
	commands["rarefaction.shared"]	= "rarefaction.shared"; 
	commands["summary.single"]		= "summary.single"; 
	commands["summary.shared"]		= "summary.shared"; 
	commands["parsimony"]			= "parsimony";
	commands["unifrac.weighted"]	= "unifrac.weighted"; 
	commands["unifrac.unweighted"]	= "unifrac.unweighted"; 
	commands["libshuff"]			= "libshuff";
	commands["tree.shared"]			= "tree.shared";
	commands["heatmap.bin"]			= "heatmap.bin";
	commands["heatmap.sim"]			= "heatmap.sim";
	commands["venn"]				= "venn";
	commands["get.group"]           = "get.group";
	commands["get.label"]           = "get.label";
	commands["get.sabund"]          = "get.sabund";
	commands["get.rabund"]          = "get.rabund";
	commands["help"]				= "help";
	commands["reverse.seqs"]		= "reverse.seqs";
	commands["trim.seqs"]			= "trim.seqs";
	commands["trim.flows"]			= "trim.flows";
	commands["list.seqs"]			= "list.seqs";
	commands["get.seqs"]			= "get.seqs";
	commands["remove.seqs"]			= "remove.seqs";
	commands["system"]				= "system";
	commands["align.check"]			= "align.check";
	commands["get.sharedseqs"]		= "get.sharedseqs";
	commands["get.otulist"]			= "get.otulist";
	commands["hcluster"]			= "hcluster"; 
	commands["phylotype"]			= "phylotype";
	commands["mgcluster"]			= "mgcluster";
	commands["pre.cluster"]			= "pre.cluster";
	commands["pcoa"]				= "pcoa";
	commands["otu.hierarchy"]		= "otu.hierarchy";
	commands["set.dir"]				= "MPIEnabled";
	commands["merge.files"]			= "merge.files";
	commands["parse.list"]			= "parse.list";
	commands["set.logfile"]			= "set.logfile";
	commands["phylo.diversity"]		= "phylo.diversity";
	commands["make.group"]			= "make.group";
	commands["chop.seqs"]			= "chop.seqs";
	commands["clearcut"]			= "clearcut";
	commands["catchall"]			= "catchall";
	commands["split.abund"]			= "split.abund";
	commands["classify.otu"]		= "classify.otu";
	commands["degap.seqs"]			= "degap.seqs";
	commands["get.relabund"]		= "get.relabund";
	commands["sffinfo"]				= "sffinfo";
	commands["normalize.shared"]	= "normalize.shared";
	commands["metastats"]			= "metastats";
	commands["split.groups"]		= "split.groups";
	commands["cluster.fragments"]	= "cluster.fragments";
	commands["get.lineage"]			= "get.lineage";
	commands["remove.lineage"]		= "remove.lineage";
	commands["fastq.info"]			= "fastq.info";
	commands["deunique.seqs"]		= "deunique.seqs";
	commands["cluster.classic"]		= "cluster.classic";
	commands["sub.sample"]			= "sub.sample";
	commands["remove.groups"]		= "remove.groups";
	commands["get.groups"]			= "get.groups";
	commands["get.otus"]			= "get.otus";
	commands["remove.otus"]			= "remove.otus";
	commands["indicator"]			= "indicator";
	commands["consensus.seqs"]		= "consensus.seqs";
	commands["corr.axes"]			= "corr.axes";
	commands["pca"]					= "pca";
	commands["nmds"]				= "nmds";
	commands["remove.rare"]			= "remove.rare";
	commands["amova"]				= "amova";
	commands["homova"]				= "homova";
	commands["mantel"]				= "mantel";
	commands["anosim"]				= "anosim";
	commands["make.fastq"]			= "make.fastq";
	commands["merge.groups"]		= "merge.groups";
	commands["get.current"]			= "MPIEnabled";
	commands["set.current"]			= "MPIEnabled";
	commands["get.commandinfo"]		= "get.commandinfo";
	commands["deunique.tree"]		= "deunique.tree";
	commands["count.seqs"]			= "count.seqs";
	commands["count.groups"]		= "count.groups";
	commands["clear.memory"]		= "clear.memory";
	commands["pairwise.seqs"]		= "MPIEnabled";
	commands["pipeline.pds"]		= "MPIEnabled";
	commands["classify.seqs"]		= "MPIEnabled"; 
	commands["dist.seqs"]			= "MPIEnabled";
	commands["filter.seqs"]			= "MPIEnabled";
	commands["align.seqs"]			= "MPIEnabled";
	commands["chimera.ccode"]		= "MPIEnabled";
	commands["chimera.check"]		= "MPIEnabled";
	commands["chimera.slayer"]		= "MPIEnabled";
	commands["chimera.uchime"]		= "chimera.uchime";
	commands["chimera.perseus"]		= "chimera.perseus";
	commands["chimera.pintail"]		= "MPIEnabled";
	commands["chimera.bellerophon"]	= "MPIEnabled";
	commands["screen.seqs"]			= "MPIEnabled";
	commands["summary.seqs"]		= "MPIEnabled";
	commands["cluster.split"]		= "MPIEnabled";
	commands["shhh.flows"]			= "MPIEnabled";
	commands["sens.spec"]			= "sens.spec";
	commands["seq.error"]			= "seq.error";
	commands["summary.tax"]			= "summary.tax";
	commands["summary.qual"]		= "summary.qual";
	commands["shhh.seqs"]			= "shhh.seqs";
	commands["otu.association"]		= "otu.association";
    commands["sort.seqs"]           = "sort.seqs";
    commands["classify.tree"]       = "classify.tree";
    commands["cooccurrence"]        = "cooccurrence";
    commands["pcr.seqs"]            = "pcr.seqs";
    commands["create.database"]     = "create.database";
    commands["make.biom"]           = "make.biom";
    commands["get.coremicrobiome"]  = "get.coremicrobiome";
    commands["list.otulabels"]      = "list.otulabels";
	commands["quit"]				= "MPIEnabled"; 

}
/***********************************************************/

/***********************************************************/
bool CommandFactory::MPIEnabled(string commandName) {
	bool mpi = false;
	it = commands.find(commandName);
	if (it != commands.end()) { 
		if (it->second == "MPIEnabled") { return true; }
	}
	return mpi;
}
/***********************************************************/

/***********************************************************/
CommandFactory::~CommandFactory(){
	_uniqueInstance = 0;
	delete command;
	delete shellcommand;
	delete pipecommand;
}
/***********************************************************/

/***********************************************************/
int CommandFactory::checkForRedirects(string optionString) {
    try {
        
        int pos = optionString.find("outputdir");
        if (pos != string::npos) { //user has set outputdir in command option string
            string outputOption = "";
            bool foundEquals = false;
            for(int i=pos;i<optionString.length();i++){
                if(optionString[i] == ',')       { break;               }		
                else if(optionString[i] == '=')  { foundEquals = true;	}
                if (foundEquals)       {   outputOption += optionString[i]; }
            }
            if(m->dirCheck(outputOption)){ 
                setOutputDirectory(outputOption); 
                m->mothurOut("Setting output directory to: " + outputOption); m->mothurOutEndLine();
            }
        }
        
        pos = optionString.find("inputdir");
        if (pos != string::npos) { //user has set inputdir in command option string
            string intputOption = "";
            bool foundEquals = false;
            for(int i=pos;i<optionString.length();i++){
                if(optionString[i] == ',')       { break;               }		
                else if(optionString[i] == '=')  { foundEquals = true;	}
                if (foundEquals)       {   intputOption += optionString[i]; }
            }
            if(m->dirCheck(intputOption)){ 
                setInputDirectory(intputOption); 
                m->mothurOut("Setting input directory to: " + intputOption); m->mothurOutEndLine();
            }
        }
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "getCommand");
		exit(1);
	}
}
/***********************************************************/

/***********************************************************/
//This function calls the appropriate command fucntions based on user input.
Command* CommandFactory::getCommand(string commandName, string optionString){
	try {
        
		delete command;   //delete the old command
        
        checkForRedirects(optionString);
        		
		//user has opted to redirect output from dir where input files are located to some other place
		if (outputDir != "") { 
			if (optionString != "") { optionString += ", outputdir=" + outputDir; }
			else { optionString += "outputdir=" + outputDir; }
		}
		
		//user has opted to redirect input from dir where mothur.exe is located to some other place
		if (inputDir != "") { 
			if (optionString != "") { optionString += ", inputdir=" + inputDir; }
			else { optionString += "inputdir=" + inputDir; }
		}
		
		if(commandName == "read.dist")					{	command = new ReadDistCommand(optionString);				}
		else if(commandName == "read.otu")				{	command = new ReadOtuCommand(optionString);					}
		else if(commandName == "read.tree")				{	command = new ReadTreeCommand(optionString);				}
		else if(commandName == "cluster")				{	command = new ClusterCommand(optionString);					}
		else if(commandName == "unique.seqs")			{	command = new DeconvoluteCommand(optionString);				}
		else if(commandName == "parsimony")				{	command = new ParsimonyCommand(optionString);				}
		else if(commandName == "help")					{	command = new HelpCommand(optionString);					}
		else if(commandName == "quit")					{	command = new QuitCommand(optionString);					}
		else if(commandName == "collect.single")		{	command = new CollectCommand(optionString);					}
		else if(commandName == "collect.shared")		{	command = new CollectSharedCommand(optionString);			}
		else if(commandName == "rarefaction.single")	{	command = new RareFactCommand(optionString);				}
		else if(commandName == "rarefaction.shared")	{	command = new RareFactSharedCommand(optionString);			}
		else if(commandName == "summary.single")		{	command = new SummaryCommand(optionString);					}
		else if(commandName == "summary.shared")		{	command = new SummarySharedCommand(optionString);			}
		else if(commandName == "unifrac.weighted")		{	command = new UnifracWeightedCommand(optionString);			}
		else if(commandName == "unifrac.unweighted")	{	command = new UnifracUnweightedCommand(optionString);		}
		else if(commandName == "get.group")             {   command = new GetgroupCommand(optionString);				}
		else if(commandName == "get.label")             {   command = new GetlabelCommand(optionString);				}
		else if(commandName == "get.sabund")            {   command = new GetSAbundCommand(optionString);				}
		else if(commandName == "get.rabund")            {   command = new GetRAbundCommand(optionString);				}
		else if(commandName == "libshuff")              {   command = new LibShuffCommand(optionString);				}
		else if(commandName == "heatmap.bin")			{   command = new HeatMapCommand(optionString);					}
		else if(commandName == "heatmap.sim")			{   command = new HeatMapSimCommand(optionString);				}
		else if(commandName == "filter.seqs")			{   command = new FilterSeqsCommand(optionString);				}
		else if(commandName == "venn")					{   command = new VennCommand(optionString);					}
		else if(commandName == "bin.seqs")				{   command = new BinSeqCommand(optionString);					}
		else if(commandName == "get.oturep")			{   command = new GetOTURepCommand(optionString);				}
		else if(commandName == "tree.shared")			{   command = new TreeGroupCommand(optionString);				}
		else if(commandName == "dist.shared")			{   command = new MatrixOutputCommand(optionString);			}
		else if(commandName == "dist.seqs")				{   command = new DistanceCommand(optionString);				}
		else if(commandName == "align.seqs")			{   command = new AlignCommand(optionString);					}
		else if(commandName == "summary.seqs")			{	command = new SeqSummaryCommand(optionString);				}
		else if(commandName == "screen.seqs")			{	command = new ScreenSeqsCommand(optionString);				}
		else if(commandName == "reverse.seqs")			{	command = new ReverseSeqsCommand(optionString);				}
		else if(commandName == "trim.seqs")				{	command = new TrimSeqsCommand(optionString);				}
		else if(commandName == "trim.flows")			{	command = new TrimFlowsCommand(optionString);				}
		else if(commandName == "shhh.flows")			{	command = new ShhherCommand(optionString);					}
		else if(commandName == "list.seqs")				{	command = new ListSeqsCommand(optionString);				}
		else if(commandName == "get.seqs")				{	command = new GetSeqsCommand(optionString);					}
		else if(commandName == "remove.seqs")			{	command = new RemoveSeqsCommand(optionString);				}
		else if(commandName == "merge.files")			{	command = new MergeFileCommand(optionString);				}
		else if(commandName == "system")				{	command = new SystemCommand(optionString);					}
		else if(commandName == "align.check")			{	command = new AlignCheckCommand(optionString);				}
		else if(commandName == "get.sharedseqs")		{	command = new GetSharedOTUCommand(optionString);			}
		else if(commandName == "get.otulist")			{	command = new GetListCountCommand(optionString);			}
		else if(commandName == "hcluster")				{	command = new HClusterCommand(optionString);				}
		else if(commandName == "classify.seqs")			{	command = new ClassifySeqsCommand(optionString);			}
		else if(commandName == "chimera.ccode")			{	command = new ChimeraCcodeCommand(optionString);			}
		else if(commandName == "chimera.check")			{	command = new ChimeraCheckCommand(optionString);			}
		else if(commandName == "chimera.slayer")		{	command = new ChimeraSlayerCommand(optionString);			}
		else if(commandName == "chimera.uchime")		{	command = new ChimeraUchimeCommand(optionString);			}
		else if(commandName == "chimera.pintail")		{	command = new ChimeraPintailCommand(optionString);			}
		else if(commandName == "chimera.bellerophon")	{	command = new ChimeraBellerophonCommand(optionString);		}
		else if(commandName == "phylotype")				{	command = new PhylotypeCommand(optionString);				}
		else if(commandName == "mgcluster")				{	command = new MGClusterCommand(optionString);				}
		else if(commandName == "pre.cluster")			{	command = new PreClusterCommand(optionString);				}
		else if(commandName == "pcoa")					{	command = new PCOACommand(optionString);					}
		else if(commandName == "pca")					{	command = new PCACommand(optionString);						}
		else if(commandName == "nmds")					{	command = new NMDSCommand(optionString);					}
		else if(commandName == "otu.hierarchy")			{	command = new OtuHierarchyCommand(optionString);			}
		else if(commandName == "set.dir")				{	command = new SetDirectoryCommand(optionString);			}
		else if(commandName == "set.logfile")			{	command = new SetLogFileCommand(optionString);				}
		else if(commandName == "parse.list")			{	command = new ParseListCommand(optionString);				}
		else if(commandName == "phylo.diversity")		{	command = new PhyloDiversityCommand(optionString);			}
		else if(commandName == "make.group")			{	command = new MakeGroupCommand(optionString);				}
		else if(commandName == "chop.seqs")				{	command = new ChopSeqsCommand(optionString);				}
		else if(commandName == "clearcut")				{	command = new ClearcutCommand(optionString);				}
		else if(commandName == "catchall")				{	command = new CatchAllCommand(optionString);				}
		else if(commandName == "split.abund")			{	command = new SplitAbundCommand(optionString);				}
		else if(commandName == "cluster.split")			{	command = new ClusterSplitCommand(optionString);			}
		else if(commandName == "classify.otu")			{	command = new ClassifyOtuCommand(optionString);				}
		else if(commandName == "degap.seqs")			{	command = new DegapSeqsCommand(optionString);				}
		else if(commandName == "get.relabund")			{	command = new GetRelAbundCommand(optionString);				}
		else if(commandName == "sens.spec")				{	command = new SensSpecCommand(optionString);				}
		else if(commandName == "seq.error")				{	command = new SeqErrorCommand(optionString);				}
		else if(commandName == "sffinfo")				{	command = new SffInfoCommand(optionString);					}
		else if(commandName == "normalize.shared")		{	command = new NormalizeSharedCommand(optionString);			}
		else if(commandName == "metastats")				{	command = new MetaStatsCommand(optionString);				}
		else if(commandName == "split.groups")			{	command = new SplitGroupCommand(optionString);				}
		else if(commandName == "cluster.fragments")		{	command = new ClusterFragmentsCommand(optionString);		}
		else if(commandName == "get.lineage")			{	command = new GetLineageCommand(optionString);				}
		else if(commandName == "remove.lineage")		{	command = new RemoveLineageCommand(optionString);			}
		else if(commandName == "get.groups")			{	command = new GetGroupsCommand(optionString);				}
		else if(commandName == "remove.groups")			{	command = new RemoveGroupsCommand(optionString);			}
		else if(commandName == "get.otus")				{	command = new GetOtusCommand(optionString);					}
		else if(commandName == "remove.otus")			{	command = new RemoveOtusCommand(optionString);				}
		else if(commandName == "fastq.info")			{	command = new ParseFastaQCommand(optionString);				}
		else if(commandName == "pipeline.pds")			{	command = new PipelineCommand(optionString);				}
		else if(commandName == "deunique.seqs")			{	command = new DeUniqueSeqsCommand(optionString);			}
		else if(commandName == "pairwise.seqs")			{	command = new PairwiseSeqsCommand(optionString);			}
		else if(commandName == "cluster.classic")		{	command = new ClusterDoturCommand(optionString);			}
		else if(commandName == "sub.sample")			{	command = new SubSampleCommand(optionString);				}
		else if(commandName == "indicator")				{	command = new IndicatorCommand(optionString);				}
		else if(commandName == "consensus.seqs")		{	command = new ConsensusSeqsCommand(optionString);			}
		else if(commandName == "corr.axes")				{	command = new CorrAxesCommand(optionString);				}
		else if(commandName == "remove.rare")			{	command = new RemoveRareCommand(optionString);				}
		else if(commandName == "merge.groups")			{	command = new MergeGroupsCommand(optionString);				}
		else if(commandName == "amova")					{	command = new AmovaCommand(optionString);					}
		else if(commandName == "homova")				{	command = new HomovaCommand(optionString);					}
		else if(commandName == "mantel")				{	command = new MantelCommand(optionString);					}
		else if(commandName == "make.fastq")			{	command = new MakeFastQCommand(optionString);				}
		else if(commandName == "get.current")			{	command = new GetCurrentCommand(optionString);				}
		else if(commandName == "set.current")			{	command = new SetCurrentCommand(optionString);				}
		else if(commandName == "anosim")				{	command = new AnosimCommand(optionString);					}
		else if(commandName == "make.shared")			{	command = new SharedCommand(optionString);					}
		else if(commandName == "get.commandinfo")		{	command = new GetCommandInfoCommand(optionString);			}
		else if(commandName == "deunique.tree")			{	command = new DeuniqueTreeCommand(optionString);			}
		else if(commandName == "count.seqs")			{	command = new CountSeqsCommand(optionString);				}
		else if(commandName == "count.groups")			{	command = new CountGroupsCommand(optionString);				}
		else if(commandName == "clear.memory")			{	command = new ClearMemoryCommand(optionString);				}
		else if(commandName == "summary.tax")			{	command = new SummaryTaxCommand(optionString);				}
		else if(commandName == "summary.qual")			{	command = new SummaryQualCommand(optionString);				}
		else if(commandName == "chimera.perseus")		{	command = new ChimeraPerseusCommand(optionString);			}
		else if(commandName == "shhh.seqs")				{	command = new ShhhSeqsCommand(optionString);				}
		else if(commandName == "otu.association")		{	command = new OTUAssociationCommand(optionString);			}
        else if(commandName == "sort.seqs")             {	command = new SortSeqsCommand(optionString);                }
        else if(commandName == "classify.tree")         {	command = new ClassifyTreeCommand(optionString);            }
        else if(commandName == "cooccurrence")          {	command = new CooccurrenceCommand(optionString);            }
        else if(commandName == "pcr.seqs")              {	command = new PcrSeqsCommand(optionString);                 }
        else if(commandName == "create.database")       {	command = new CreateDatabaseCommand(optionString);          }
        else if(commandName == "make.biom")             {	command = new MakeBiomCommand(optionString);                }
        else if(commandName == "get.coremicrobiome")    {	command = new GetCoreMicroBiomeCommand(optionString);       }
        else if(commandName == "list.otulabels")        {	command = new ListOtuLabelsCommand(optionString);           }
		else											{	command = new NoCommand(optionString);						}

		return command;
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "getCommand");
		exit(1);
	}
}
/***********************************************************/

/***********************************************************/
//This function calls the appropriate command fucntions based on user input.
Command* CommandFactory::getCommand(string commandName, string optionString, string mode){
	try {
		delete pipecommand;   //delete the old command
		
        checkForRedirects(optionString);
        
		//user has opted to redirect output from dir where input files are located to some other place
		if (outputDir != "") { 
			if (optionString != "") { optionString += ", outputdir=" + outputDir; }
			else { optionString += "outputdir=" + outputDir; }
		}
		
		//user has opted to redirect input from dir where mothur.exe is located to some other place
		if (inputDir != "") { 
			if (optionString != "") { optionString += ", inputdir=" + inputDir; }
			else { optionString += "inputdir=" + inputDir; }
		}
		
		if(commandName == "read.dist")					{	pipecommand = new ReadDistCommand(optionString);				}
		else if(commandName == "read.otu")				{	pipecommand = new ReadOtuCommand(optionString);					}
		else if(commandName == "read.tree")				{	pipecommand = new ReadTreeCommand(optionString);				}
		else if(commandName == "cluster")				{	pipecommand = new ClusterCommand(optionString);					}
		else if(commandName == "unique.seqs")			{	pipecommand = new DeconvoluteCommand(optionString);				}
		else if(commandName == "parsimony")				{	pipecommand = new ParsimonyCommand(optionString);				}
		else if(commandName == "help")					{	pipecommand = new HelpCommand(optionString);					}
		else if(commandName == "quit")					{	pipecommand = new QuitCommand(optionString);					}
		else if(commandName == "collect.single")		{	pipecommand = new CollectCommand(optionString);					}
		else if(commandName == "collect.shared")		{	pipecommand = new CollectSharedCommand(optionString);			}
		else if(commandName == "rarefaction.single")	{	pipecommand = new RareFactCommand(optionString);				}
		else if(commandName == "rarefaction.shared")	{	pipecommand = new RareFactSharedCommand(optionString);			}
		else if(commandName == "summary.single")		{	pipecommand = new SummaryCommand(optionString);					}
		else if(commandName == "summary.shared")		{	pipecommand = new SummarySharedCommand(optionString);			}
		else if(commandName == "unifrac.weighted")		{	pipecommand = new UnifracWeightedCommand(optionString);			}
		else if(commandName == "unifrac.unweighted")	{	pipecommand = new UnifracUnweightedCommand(optionString);		}
		else if(commandName == "get.group")             {   pipecommand = new GetgroupCommand(optionString);				}
		else if(commandName == "get.label")             {   pipecommand = new GetlabelCommand(optionString);				}
		else if(commandName == "get.sabund")            {   pipecommand = new GetSAbundCommand(optionString);				}
		else if(commandName == "get.rabund")            {   pipecommand = new GetRAbundCommand(optionString);				}
		else if(commandName == "libshuff")              {   pipecommand = new LibShuffCommand(optionString);				}
		else if(commandName == "heatmap.bin")			{   pipecommand = new HeatMapCommand(optionString);					}
		else if(commandName == "heatmap.sim")			{   pipecommand = new HeatMapSimCommand(optionString);				}
		else if(commandName == "filter.seqs")			{   pipecommand = new FilterSeqsCommand(optionString);				}
		else if(commandName == "venn")					{   pipecommand = new VennCommand(optionString);					}
		else if(commandName == "bin.seqs")				{   pipecommand = new BinSeqCommand(optionString);					}
		else if(commandName == "get.oturep")			{   pipecommand = new GetOTURepCommand(optionString);				}
		else if(commandName == "tree.shared")			{   pipecommand = new TreeGroupCommand(optionString);				}
		else if(commandName == "dist.shared")			{   pipecommand = new MatrixOutputCommand(optionString);			}
		else if(commandName == "dist.seqs")				{   pipecommand = new DistanceCommand(optionString);				}
		else if(commandName == "align.seqs")			{   pipecommand = new AlignCommand(optionString);					}
		else if(commandName == "summary.seqs")			{	pipecommand = new SeqSummaryCommand(optionString);				}
		else if(commandName == "screen.seqs")			{	pipecommand = new ScreenSeqsCommand(optionString);				}
		else if(commandName == "reverse.seqs")			{	pipecommand = new ReverseSeqsCommand(optionString);				}
		else if(commandName == "trim.seqs")				{	pipecommand = new TrimSeqsCommand(optionString);				}
		else if(commandName == "trim.flows")			{	pipecommand = new TrimFlowsCommand(optionString);				}
		else if(commandName == "shhh.flows")			{	pipecommand = new ShhherCommand(optionString);					}
		else if(commandName == "list.seqs")				{	pipecommand = new ListSeqsCommand(optionString);				}
		else if(commandName == "get.seqs")				{	pipecommand = new GetSeqsCommand(optionString);					}
		else if(commandName == "remove.seqs")			{	pipecommand = new RemoveSeqsCommand(optionString);				}
		else if(commandName == "merge.files")			{	pipecommand = new MergeFileCommand(optionString);				}
		else if(commandName == "system")				{	pipecommand = new SystemCommand(optionString);					}
		else if(commandName == "align.check")			{	pipecommand = new AlignCheckCommand(optionString);				}
		else if(commandName == "get.sharedseqs")		{	pipecommand = new GetSharedOTUCommand(optionString);			}
		else if(commandName == "get.otulist")			{	pipecommand = new GetListCountCommand(optionString);			}
		else if(commandName == "hcluster")				{	pipecommand = new HClusterCommand(optionString);				}
		else if(commandName == "classify.seqs")			{	pipecommand = new ClassifySeqsCommand(optionString);			}
		else if(commandName == "chimera.ccode")			{	pipecommand = new ChimeraCcodeCommand(optionString);			}
		else if(commandName == "chimera.check")			{	pipecommand = new ChimeraCheckCommand(optionString);			}
		else if(commandName == "chimera.uchime")		{	pipecommand = new ChimeraUchimeCommand(optionString);			}
		else if(commandName == "chimera.slayer")		{	pipecommand = new ChimeraSlayerCommand(optionString);			}
		else if(commandName == "chimera.pintail")		{	pipecommand = new ChimeraPintailCommand(optionString);			}
		else if(commandName == "chimera.bellerophon")	{	pipecommand = new ChimeraBellerophonCommand(optionString);		}
		else if(commandName == "phylotype")				{	pipecommand = new PhylotypeCommand(optionString);				}
		else if(commandName == "mgcluster")				{	pipecommand = new MGClusterCommand(optionString);				}
		else if(commandName == "pre.cluster")			{	pipecommand = new PreClusterCommand(optionString);				}
		else if(commandName == "pcoa")					{	pipecommand = new PCOACommand(optionString);					}
		else if(commandName == "pca")					{	pipecommand = new PCACommand(optionString);						}
		else if(commandName == "nmds")					{	pipecommand = new NMDSCommand(optionString);					}
		else if(commandName == "otu.hierarchy")			{	pipecommand = new OtuHierarchyCommand(optionString);			}
		else if(commandName == "set.dir")				{	pipecommand = new SetDirectoryCommand(optionString);			}
		else if(commandName == "set.logfile")			{	pipecommand = new SetLogFileCommand(optionString);				}
		else if(commandName == "parse.list")			{	pipecommand = new ParseListCommand(optionString);				}
		else if(commandName == "phylo.diversity")		{	pipecommand = new PhyloDiversityCommand(optionString);			}
		else if(commandName == "make.group")			{	pipecommand = new MakeGroupCommand(optionString);				}
		else if(commandName == "chop.seqs")				{	pipecommand = new ChopSeqsCommand(optionString);				}
		else if(commandName == "clearcut")				{	pipecommand = new ClearcutCommand(optionString);				}
		else if(commandName == "catchall")				{	pipecommand = new CatchAllCommand(optionString);				}
		else if(commandName == "split.abund")			{	pipecommand = new SplitAbundCommand(optionString);				}
		else if(commandName == "cluster.split")			{	pipecommand = new ClusterSplitCommand(optionString);			}
		else if(commandName == "classify.otu")			{	pipecommand = new ClassifyOtuCommand(optionString);				}
		else if(commandName == "degap.seqs")			{	pipecommand = new DegapSeqsCommand(optionString);				}
		else if(commandName == "get.relabund")			{	pipecommand = new GetRelAbundCommand(optionString);				}
		else if(commandName == "sens.spec")				{	pipecommand = new SensSpecCommand(optionString);				}
		else if(commandName == "seq.error")				{	pipecommand = new SeqErrorCommand(optionString);				}
		else if(commandName == "sffinfo")				{	pipecommand = new SffInfoCommand(optionString);					}
		else if(commandName == "normalize.shared")		{	pipecommand = new NormalizeSharedCommand(optionString);			}
		else if(commandName == "metastats")				{	pipecommand = new MetaStatsCommand(optionString);				}
		else if(commandName == "split.groups")			{	pipecommand = new SplitGroupCommand(optionString);				}
		else if(commandName == "cluster.fragments")		{	pipecommand = new ClusterFragmentsCommand(optionString);		}
		else if(commandName == "get.lineage")			{	pipecommand = new GetLineageCommand(optionString);				}
		else if(commandName == "get.groups")			{	pipecommand = new GetGroupsCommand(optionString);				}
		else if(commandName == "remove.lineage")		{	pipecommand = new RemoveLineageCommand(optionString);			}
		else if(commandName == "remove.groups")			{	pipecommand = new RemoveGroupsCommand(optionString);			}
		else if(commandName == "get.otus")				{	pipecommand = new GetOtusCommand(optionString);					}
		else if(commandName == "remove.otus")			{	pipecommand = new RemoveOtusCommand(optionString);				}
		else if(commandName == "fastq.info")			{	pipecommand = new ParseFastaQCommand(optionString);				}
		else if(commandName == "deunique.seqs")			{	pipecommand = new DeUniqueSeqsCommand(optionString);			}
		else if(commandName == "pairwise.seqs")			{	pipecommand = new PairwiseSeqsCommand(optionString);			}
		else if(commandName == "cluster.classic")		{	pipecommand = new ClusterDoturCommand(optionString);			}
		else if(commandName == "sub.sample")			{	pipecommand = new SubSampleCommand(optionString);				}
		else if(commandName == "indicator")				{	pipecommand = new IndicatorCommand(optionString);				}
		else if(commandName == "consensus.seqs")		{	pipecommand = new ConsensusSeqsCommand(optionString);			}
		else if(commandName == "corr.axes")				{	pipecommand = new CorrAxesCommand(optionString);				}
		else if(commandName == "remove.rare")			{	pipecommand = new RemoveRareCommand(optionString);				}
		else if(commandName == "merge.groups")			{	pipecommand = new MergeGroupsCommand(optionString);				}
		else if(commandName == "amova")					{	pipecommand = new AmovaCommand(optionString);					}
		else if(commandName == "homova")				{	pipecommand = new HomovaCommand(optionString);					}
		else if(commandName == "mantel")				{	pipecommand = new MantelCommand(optionString);					}
		else if(commandName == "anosim")				{	pipecommand = new AnosimCommand(optionString);					}
		else if(commandName == "make.fastq")			{	pipecommand = new MakeFastQCommand(optionString);				}
		else if(commandName == "get.current")			{	pipecommand = new GetCurrentCommand(optionString);				}
		else if(commandName == "set.current")			{	pipecommand = new SetCurrentCommand(optionString);				}
		else if(commandName == "make.shared")			{	pipecommand = new SharedCommand(optionString);					}
		else if(commandName == "get.commandinfo")		{	pipecommand = new GetCommandInfoCommand(optionString);			}
		else if(commandName == "deunique.tree")			{	pipecommand = new DeuniqueTreeCommand(optionString);			}
		else if(commandName == "count.seqs")			{	pipecommand = new CountSeqsCommand(optionString);				}
		else if(commandName == "count.groups")			{	pipecommand = new CountGroupsCommand(optionString);				}
		else if(commandName == "clear.memory")			{	pipecommand = new ClearMemoryCommand(optionString);				}
		else if(commandName == "summary.tax")			{	pipecommand = new SummaryTaxCommand(optionString);				}
		else if(commandName == "summary.qual")			{	pipecommand = new SummaryQualCommand(optionString);				}
		else if(commandName == "chimera.perseus")		{	pipecommand = new ChimeraPerseusCommand(optionString);			}
		else if(commandName == "shhh.seqs")				{	pipecommand = new ShhhSeqsCommand(optionString);				}
		else if(commandName == "otu.association")		{	pipecommand = new OTUAssociationCommand(optionString);			}
        else if(commandName == "sort.seqs")             {	pipecommand = new SortSeqsCommand(optionString);                }
        else if(commandName == "classify.tree")         {	pipecommand = new ClassifyTreeCommand(optionString);            }
        else if(commandName == "cooccurrence")          {	pipecommand = new CooccurrenceCommand(optionString);            }
        else if(commandName == "pcr.seqs")              {	pipecommand = new PcrSeqsCommand(optionString);                 }
        else if(commandName == "create.database")       {	pipecommand = new CreateDatabaseCommand(optionString);          }
        else if(commandName == "make.biom")             {	pipecommand = new MakeBiomCommand(optionString);                }
        else if(commandName == "get.coremicrobiome")    {	pipecommand = new GetCoreMicroBiomeCommand(optionString);       }
        else if(commandName == "list.otulabels")        {	pipecommand = new ListOtuLabelsCommand(optionString);           }
		else											{	pipecommand = new NoCommand(optionString);						}

		return pipecommand;
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "getCommand");
		exit(1);
	}
}
/***********************************************************/

/***********************************************************/
//This function calls the appropriate command fucntions based on user input, this is used by the pipeline command to check a users piepline for errors before running
Command* CommandFactory::getCommand(string commandName){
	try {
		delete shellcommand;   //delete the old command
		
		if(commandName == "read.dist")					{	shellcommand = new ReadDistCommand();				}
		else if(commandName == "read.otu")				{	shellcommand = new ReadOtuCommand();				}
		else if(commandName == "read.tree")				{	shellcommand = new ReadTreeCommand();				}
		else if(commandName == "cluster")				{	shellcommand = new ClusterCommand();				}
		else if(commandName == "unique.seqs")			{	shellcommand = new DeconvoluteCommand();			}
		else if(commandName == "parsimony")				{	shellcommand = new ParsimonyCommand();				}
		else if(commandName == "help")					{	shellcommand = new HelpCommand();					}
		else if(commandName == "quit")					{	shellcommand = new QuitCommand();					}
		else if(commandName == "collect.single")		{	shellcommand = new CollectCommand();				}
		else if(commandName == "collect.shared")		{	shellcommand = new CollectSharedCommand();			}
		else if(commandName == "rarefaction.single")	{	shellcommand = new RareFactCommand();				}
		else if(commandName == "rarefaction.shared")	{	shellcommand = new RareFactSharedCommand();			}
		else if(commandName == "summary.single")		{	shellcommand = new SummaryCommand();				}
		else if(commandName == "summary.shared")		{	shellcommand = new SummarySharedCommand();			}
		else if(commandName == "unifrac.weighted")		{	shellcommand = new UnifracWeightedCommand();		}
		else if(commandName == "unifrac.unweighted")	{	shellcommand = new UnifracUnweightedCommand();		}
		else if(commandName == "get.group")             {   shellcommand = new GetgroupCommand();				}
		else if(commandName == "get.label")             {   shellcommand = new GetlabelCommand();				}
		else if(commandName == "get.sabund")            {   shellcommand = new GetSAbundCommand();				}
		else if(commandName == "get.rabund")            {   shellcommand = new GetRAbundCommand();				}
		else if(commandName == "libshuff")              {   shellcommand = new LibShuffCommand();				}
		else if(commandName == "heatmap.bin")			{   shellcommand = new HeatMapCommand();				}
		else if(commandName == "heatmap.sim")			{   shellcommand = new HeatMapSimCommand();				}
		else if(commandName == "filter.seqs")			{   shellcommand = new FilterSeqsCommand();				}
		else if(commandName == "venn")					{   shellcommand = new VennCommand();					}
		else if(commandName == "bin.seqs")				{   shellcommand = new BinSeqCommand();					}
		else if(commandName == "get.oturep")			{   shellcommand = new GetOTURepCommand();				}
		else if(commandName == "tree.shared")			{   shellcommand = new TreeGroupCommand();				}
		else if(commandName == "dist.shared")			{   shellcommand = new MatrixOutputCommand();			}
		else if(commandName == "dist.seqs")				{   shellcommand = new DistanceCommand();				}
		else if(commandName == "align.seqs")			{   shellcommand = new AlignCommand();					}
		else if(commandName == "summary.seqs")			{	shellcommand = new SeqSummaryCommand();				}
		else if(commandName == "screen.seqs")			{	shellcommand = new ScreenSeqsCommand();				}
		else if(commandName == "reverse.seqs")			{	shellcommand = new ReverseSeqsCommand();			}
		else if(commandName == "trim.seqs")				{	shellcommand = new TrimSeqsCommand();				}
		else if(commandName == "trim.flows")			{	shellcommand = new TrimFlowsCommand();				}
		else if(commandName == "shhh.flows")			{	shellcommand = new ShhherCommand();					}
		else if(commandName == "list.seqs")				{	shellcommand = new ListSeqsCommand();				}
		else if(commandName == "get.seqs")				{	shellcommand = new GetSeqsCommand();				}
		else if(commandName == "remove.seqs")			{	shellcommand = new RemoveSeqsCommand();				}
		else if(commandName == "merge.files")			{	shellcommand = new MergeFileCommand();				}
		else if(commandName == "system")				{	shellcommand = new SystemCommand();					}
		else if(commandName == "align.check")			{	shellcommand = new AlignCheckCommand();				}
		else if(commandName == "get.sharedseqs")		{	shellcommand = new GetSharedOTUCommand();			}
		else if(commandName == "get.otulist")			{	shellcommand = new GetListCountCommand();			}
		else if(commandName == "hcluster")				{	shellcommand = new HClusterCommand();				}
		else if(commandName == "classify.seqs")			{	shellcommand = new ClassifySeqsCommand();			}
		else if(commandName == "chimera.ccode")			{	shellcommand = new ChimeraCcodeCommand();			}
		else if(commandName == "chimera.check")			{	shellcommand = new ChimeraCheckCommand();			}
		else if(commandName == "chimera.slayer")		{	shellcommand = new ChimeraSlayerCommand();			}
		else if(commandName == "chimera.uchime")		{	shellcommand = new ChimeraUchimeCommand();			}
		else if(commandName == "chimera.pintail")		{	shellcommand = new ChimeraPintailCommand();			}
		else if(commandName == "chimera.bellerophon")	{	shellcommand = new ChimeraBellerophonCommand();		}
		else if(commandName == "phylotype")				{	shellcommand = new PhylotypeCommand();				}
		else if(commandName == "mgcluster")				{	shellcommand = new MGClusterCommand();				}
		else if(commandName == "pre.cluster")			{	shellcommand = new PreClusterCommand();				}
		else if(commandName == "pcoa")					{	shellcommand = new PCOACommand();					}
		else if(commandName == "pca")					{	shellcommand = new PCACommand();					}
		else if(commandName == "nmds")					{	shellcommand = new NMDSCommand();					}
		else if(commandName == "otu.hierarchy")			{	shellcommand = new OtuHierarchyCommand();			}
		else if(commandName == "set.dir")				{	shellcommand = new SetDirectoryCommand();			}
		else if(commandName == "set.logfile")			{	shellcommand = new SetLogFileCommand();				}
		else if(commandName == "parse.list")			{	shellcommand = new ParseListCommand();				}
		else if(commandName == "phylo.diversity")		{	shellcommand = new PhyloDiversityCommand();			}
		else if(commandName == "make.group")			{	shellcommand = new MakeGroupCommand();				}
		else if(commandName == "chop.seqs")				{	shellcommand = new ChopSeqsCommand();				}
		else if(commandName == "clearcut")				{	shellcommand = new ClearcutCommand();				}
		else if(commandName == "catchall")				{	shellcommand = new CatchAllCommand();				}
		else if(commandName == "split.abund")			{	shellcommand = new SplitAbundCommand();				}
		else if(commandName == "cluster.split")			{	shellcommand = new ClusterSplitCommand();			}
		else if(commandName == "classify.otu")			{	shellcommand = new ClassifyOtuCommand();			}
		else if(commandName == "degap.seqs")			{	shellcommand = new DegapSeqsCommand();				}
		else if(commandName == "get.relabund")			{	shellcommand = new GetRelAbundCommand();			}
		else if(commandName == "sens.spec")				{	shellcommand = new SensSpecCommand();				}
		else if(commandName == "seq.error")				{	shellcommand = new SeqErrorCommand();				}
		else if(commandName == "sffinfo")				{	shellcommand = new SffInfoCommand();				}
		else if(commandName == "normalize.shared")		{	shellcommand = new NormalizeSharedCommand();		}
		else if(commandName == "metastats")				{	shellcommand = new MetaStatsCommand();				}
		else if(commandName == "split.groups")			{	shellcommand = new SplitGroupCommand();				}
		else if(commandName == "cluster.fragments")		{	shellcommand = new ClusterFragmentsCommand();		}
		else if(commandName == "get.lineage")			{	shellcommand = new GetLineageCommand();				}
		else if(commandName == "remove.lineage")		{	shellcommand = new RemoveLineageCommand();			}
		else if(commandName == "get.groups")			{	shellcommand = new GetGroupsCommand();				}
		else if(commandName == "remove.groups")			{	shellcommand = new RemoveGroupsCommand();			}
		else if(commandName == "get.otus")				{	shellcommand = new GetOtusCommand();				}
		else if(commandName == "remove.otus")			{	shellcommand = new RemoveOtusCommand();				}
		else if(commandName == "fastq.info")			{	shellcommand = new ParseFastaQCommand();			}
		else if(commandName == "deunique.seqs")			{	shellcommand = new DeUniqueSeqsCommand();			}
		else if(commandName == "pairwise.seqs")			{	shellcommand = new PairwiseSeqsCommand();			}
		else if(commandName == "cluster.classic")		{	shellcommand = new ClusterDoturCommand();			}
		else if(commandName == "sub.sample")			{	shellcommand = new SubSampleCommand();				}
		else if(commandName == "indicator")				{	shellcommand = new IndicatorCommand();				}
		else if(commandName == "consensus.seqs")		{	shellcommand = new ConsensusSeqsCommand();			}
		else if(commandName == "corr.axes")				{	shellcommand = new CorrAxesCommand();				}
		else if(commandName == "remove.rare")			{	shellcommand = new RemoveRareCommand();				}
		else if(commandName == "merge.groups")			{	shellcommand = new MergeGroupsCommand();			}
		else if(commandName == "amova")					{	shellcommand = new AmovaCommand();					}
		else if(commandName == "homova")				{	shellcommand = new HomovaCommand();					}
		else if(commandName == "mantel")				{	shellcommand = new MantelCommand();					}
		else if(commandName == "anosim")				{	shellcommand = new AnosimCommand();					}
		else if(commandName == "make.fastq")			{	shellcommand = new MakeFastQCommand();				}
		else if(commandName == "get.current")			{	shellcommand = new GetCurrentCommand();				}
		else if(commandName == "set.current")			{	shellcommand = new SetCurrentCommand();				}
		else if(commandName == "make.shared")			{	shellcommand = new SharedCommand();					}
		else if(commandName == "get.commandinfo")		{	shellcommand = new GetCommandInfoCommand();			}
		else if(commandName == "deunique.tree")			{	shellcommand = new DeuniqueTreeCommand();			}
		else if(commandName == "count.seqs")			{	shellcommand = new CountSeqsCommand();				}
		else if(commandName == "count.groups")			{	shellcommand = new CountGroupsCommand();			}
		else if(commandName == "clear.memory")			{	shellcommand = new ClearMemoryCommand();			}
		else if(commandName == "summary.tax")			{	shellcommand = new SummaryTaxCommand();				}
		else if(commandName == "summary.qual")			{	shellcommand = new SummaryQualCommand();			}
		else if(commandName == "chimera.perseus")		{	shellcommand = new ChimeraPerseusCommand();			}
		else if(commandName == "shhh.seqs")				{	shellcommand = new ShhhSeqsCommand();				}
		else if(commandName == "otu.association")		{	shellcommand = new OTUAssociationCommand();			}
        else if(commandName == "sort.seqs")             {	shellcommand = new SortSeqsCommand();               }
        else if(commandName == "classify.tree")         {	shellcommand = new ClassifyTreeCommand();           }
        else if(commandName == "cooccurrence")          {	shellcommand = new CooccurrenceCommand();           }
        else if(commandName == "pcr.seqs")              {	shellcommand = new PcrSeqsCommand();                }
        else if(commandName == "create.database")       {	shellcommand = new CreateDatabaseCommand();         }
        else if(commandName == "make.biom")             {	shellcommand = new MakeBiomCommand();               }
        else if(commandName == "get.coremicrobiome")    {	shellcommand = new GetCoreMicroBiomeCommand();      }
        else if(commandName == "list.otulabels")        {	shellcommand = new ListOtuLabelsCommand();          }
		else											{	shellcommand = new NoCommand();						}

		return shellcommand;
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "getCommand");
		exit(1);
	}
}
/***********************************************************
//This function is used to interrupt a command
Command* CommandFactory::getCommand(){
	try {
		delete command;   //delete the old command

		string s = "";
	    command = new NoCommand(s);
	
		return command;
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "getCommand");
		exit(1);
	}
}
***********************************************************************/
bool CommandFactory::isValidCommand(string command) {
	try {	
	
		//is the command in the map
		if ((commands.find(command)) != (commands.end())) {
			return true;
		}else{
			m->mothurOut(command + " is not a valid command in Mothur.  Valid commands are ");
			for (it = commands.begin(); it != commands.end(); it++) {
				m->mothurOut(it->first + ", ");
			}
			m->mothurOutEndLine();
			return false;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "isValidCommand");
		exit(1);
	}
}
/***********************************************************************/
bool CommandFactory::isValidCommand(string command, string noError) {
	try {	
	
		//is the command in the map
		if ((commands.find(command)) != (commands.end())) {
			return true;
		}else{
			return false;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "isValidCommand");
		exit(1);
	}
}
/***********************************************************************/
void CommandFactory::printCommands(ostream& out) {
	try {	
		out << "Valid commands are: ";
		for (it = commands.begin(); it != commands.end(); it++) {
			out << it->first << ",";
		}
		out << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "CommandFactory", "printCommands");
		exit(1);
	}
}
/***********************************************************************/




