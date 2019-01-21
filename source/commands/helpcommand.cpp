/*
 *  helpcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "helpcommand.h"
#include "command.hpp"
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
#include "deuniqueseqscommand.h"
#include "pairwiseseqscommand.h"
#include "clusterdoturcommand.h"
#include "subsamplecommand.h"
#include "removegroupscommand.h"
#include "getgroupscommand.h"
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
#include "getotulabelscommand.h"
#include "removeotulabelscommand.h"
#include "makecontigscommand.h"
#include "sffmultiplecommand.h"
#include "classifysvmsharedcommand.h"
#include "classifyrfsharedcommand.h"
#include "filtersharedcommand.h"
#include "primerdesigncommand.h"
#include "getdistscommand.h"
#include "removedistscommand.h"
#include "mergetaxsummarycommand.h"
#include "getmetacommunitycommand.h"
#include "sparcccommand.h"
#include "makelookupcommand.h"
#include "renameseqscommand.h"
#include "makelefsecommand.h"
#include "lefsecommand.h"
#include "kruskalwalliscommand.h"
#include "sracommand.h"
#include "mergesfffilecommand.h"
#include "getmimarkspackagecommand.h"
#include "mimarksattributescommand.h"
#include "setseedcommand.h"
#include "makefilecommand.h"
#include "biominfocommand.h"
#include "renamefilecommand.h"
#include "chimeravsearchcommand.h"
#include "mergecountcommand.hpp"


//**********************************************************************************************************************

HelpCommand::HelpCommand(string option)  {	
	validCommands = CommandFactory::getInstance();
    
    abort = false; calledHelp = false;
    
    //allow user to run help
    if(option == "help") { help(); abort = true; calledHelp = true; }
    else if(option == "citation") { citation(); abort = true; calledHelp = true;}
    
    commandName = option;
    
}
//**********************************************************************************************************************
int HelpCommand::execute(){
	try {
        if (commandName != "") {
            if (validCommands->isValidCommand(commandName)) {
                Command* command;
                string optionString = "help";
                
                if(commandName == "cluster")                    {	command = new ClusterCommand(optionString);					}
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
                else if(commandName == "classify.seqs")			{	command = new ClassifySeqsCommand(optionString);			}
                else if(commandName == "chimera.ccode")			{	command = new ChimeraCcodeCommand(optionString);			}
                else if(commandName == "chimera.check")			{	command = new ChimeraCheckCommand(optionString);			}
                else if(commandName == "chimera.slayer")		{	command = new ChimeraSlayerCommand(optionString);			}
                else if(commandName == "chimera.uchime")		{	command = new ChimeraUchimeCommand(optionString);			}
                else if(commandName == "chimera.pintail")		{	command = new ChimeraPintailCommand(optionString);			}
                else if(commandName == "chimera.bellerophon")	{	command = new ChimeraBellerophonCommand(optionString);		}
                else if(commandName == "chimera.vsearch")       {	command = new ChimeraVsearchCommand(optionString);          }
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
                else if((commandName == "get.otus")	|| (commandName == "get.otulabels"))			{	command = new GetOtuLabelsCommand(optionString);			}
                else if((commandName == "remove.otus") || (commandName == "remove.otulabels"))			{	command = new RemoveOtuLabelsCommand(optionString);			}
                else if((commandName == "list.otus")	||(commandName == "list.otulabels"))        {	command = new ListOtuLabelsCommand(optionString);           }
                else if(commandName == "fastq.info")			{	command = new ParseFastaQCommand(optionString);				}
                else if(commandName == "deunique.seqs")			{	command = new DeUniqueSeqsCommand(optionString);			}
                else if(commandName == "pairwise.seqs")			{	command = new PairwiseSeqsCommand(optionString);			}
                else if(commandName == "cluster.classic")		{	command = new ClusterDoturCommand(optionString);			}
                else if(commandName == "sub.sample")			{	command = new SubSampleCommand(optionString);				}
                else if(commandName == "indicator")				{	command = new IndicatorCommand(optionString);				}
                else if(commandName == "consensus.seqs")		{	command = new ConsensusSeqsCommand(optionString);			}
                else if(commandName == "corr.axes")				{	command = new CorrAxesCommand(optionString);				}
                else if(commandName == "remove.rare")			{	command = new RemoveRareCommand(optionString);				}
                else if(commandName == "merge.groups")			{	command = new MergeGroupsCommand(optionString);				}
                else if(commandName == "merge.count")			{	command = new MergeCountCommand(optionString);				}
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
                else if((commandName == "count.seqs") || (commandName == "make.table"))			{	command = new CountSeqsCommand(optionString);				}
                else if(commandName == "count.groups")			{	command = new CountGroupsCommand(optionString);				}
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
                else if(commandName == "make.contigs")          {	command = new MakeContigsCommand(optionString);             }
                else if(commandName == "sff.multiple")          {	command = new SffMultipleCommand(optionString);             }
                else if(commandName == "classify.svm")          {   command = new ClassifySvmSharedCommand(optionString);       }
                else if(commandName == "classify.rf")           {	command = new ClassifyRFSharedCommand(optionString);          }
                else if(commandName == "filter.shared")         {	command = new FilterSharedCommand(optionString);            }
                else if(commandName == "primer.design")         {	command = new PrimerDesignCommand(optionString);            }
                else if(commandName == "get.dists")             {	command = new GetDistsCommand(optionString);                }
                else if(commandName == "remove.dists")          {	command = new RemoveDistsCommand(optionString);             }
                else if(commandName == "merge.taxsummary")      {	command = new MergeTaxSummaryCommand(optionString);         }
                else if(commandName == "get.communitytype")     {	command = new GetMetaCommunityCommand(optionString);        }
                else if(commandName == "sparcc")                {	command = new SparccCommand(optionString);                  }
                else if(commandName == "make.lookup")			{	command = new MakeLookupCommand(optionString);				}
                else if(commandName == "rename.seqs")			{	command = new RenameSeqsCommand(optionString);				}
                else if(commandName == "make.lefse")			{	command = new MakeLefseCommand(optionString);				}
                else if(commandName == "lefse")                 {	command = new LefseCommand(optionString);                   }
                else if(commandName == "kruskal.wallis")        {	command = new KruskalWallisCommand(optionString);           }
                else if(commandName == "make.sra")              {	command = new SRACommand(optionString);                     }
                else if(commandName == "merge.sfffiles")        {	command = new MergeSfffilesCommand(optionString);           }
                else if(commandName == "get.mimarkspackage")    {	command = new GetMIMarksPackageCommand(optionString);       }
                else if(commandName == "mimarks.attributes")    {	command = new MimarksAttributesCommand(optionString);       }
                else if(commandName == "set.seed")              {	command = new SetSeedCommand(optionString);                 }
                else if(commandName == "make.file")             {	command = new MakeFileCommand(optionString);                }
                else if(commandName == "biom.info")             {	command = new BiomInfoCommand(optionString);                }
                else if(commandName == "rename.file")           {	command = new RenameFileCommand(optionString);              }
                else											{	command = new NoCommand(optionString);						}
                
                command->execute();
                delete command;
                
            }else {
                m->mothurOut("[ERROR]: " + commandName + " is not a valid command."); m->mothurOutEndLine();
                validCommands->printCommands(cout);
            }
        }else {
            validCommands->printCommands(cout);
            m->mothurOut("For more information about a specific command type 'commandName(help)' i.e. 'cluster(help)'"); m->mothurOutEndLine();
        }
		
        m->mothurOutEndLine(); m->mothurOut("For further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki, or contact Pat Schloss at mothur.bugs@gmail.com.\n");
	
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HelpCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************/
