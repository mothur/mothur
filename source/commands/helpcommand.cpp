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
#include "uniqueseqscommand.h"
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
#include "treesharedcommand.h"
#include "distancecommand.h"
#include "aligncommand.h"
#include "distsharedcommand.h"
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
#include "aligncheckcommand.h"
#include "getsharedotucommand.h"
#include "getlistcountcommand.h"
#include "classifyseqscommand.h"
#include "phylotypecommand.h"
#include "mgclustercommand.h"
#include "preclustercommand.h"
#include "pcoacommand.h"
#include "otuhierarchycommand.h"
#include "setdircommand.h"
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
#include "makesharedcommand.h"
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
#include "listotuscommand.h"
#include "getotuscommand.h"
#include "removeotuscommand.h"
#include "makecontigscommand.h"
#include "sffmultiplecommand.h"
#include "classifysvmsharedcommand.h"
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

HelpCommand::HelpCommand(string option) : Command()  {	
	validCommands = CommandFactory::getInstance();
    
    abort = false; calledHelp = false;
    
    //allow user to run help
    if(option == "help") { help(); abort = true; calledHelp = true; }
    else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
    
    commandName = option;
    
}
//**********************************************************************************************************************
string HelpCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
    
        string question = "How do I cite mothur?"; questions.push_back(question);
        string qanswer = "\tSchloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.\n"; qanswers.push_back(qanswer);
        
        question = "Do you have an example analysis?"; questions.push_back(question);
        qanswer = "\tYes, https://mothur.org/wiki/454_SOP and https://mothur.org/wiki/MiSeq_SOP highlight some of the things you can do with mothur.\n"; qanswers.push_back(qanswer);
        
        question = "Do you offer workshops?"; questions.push_back(question);
        qanswer = "\tYes! Please see our https://mothur.org/wiki/Workshops page for more information.\n"; qanswers.push_back(qanswer);
        
        question = "What are mothur's file types?"; questions.push_back(question);
        qanswer = "\tMothur uses and creates many file types. Including fasta, name, group, design, count, list, rabund, sabund, shared, relabund, oligos, taxonomy, constaxonomy, phylip, column, flow, qfile, file, biom and tree. You can find out more about these formats here: https://www.mothur.org/wiki/File_Types.\n"; qanswers.push_back(qanswer);
        
        question = "Is there a list of all of mothur's commands?"; questions.push_back(question);
        qanswer = "\tYes! You can find it here, http://www.mothur.org/wiki/Category:Commands.\n"; qanswers.push_back(qanswer);
        
        question = "Why does the cutoff change when I cluster with average neighbor?"; questions.push_back(question);
        qanswer = "\tThis is a product of using the average neighbor algorithm with a sparse distance matrix. When you run cluster, the algorithm looks for pairs of sequences to merge in the rows and columns that are getting merged together. Let's say you set the cutoff to 0.05. If one cell has a distance of 0.03 and the cell it is getting merged with has a distance above 0.05 then the cutoff is reset to 0.03, because it's not possible to merge at a higher level and keep all the data. All of the sequences are still there from multiple phyla. Incidentally, although we always see this, it is a bigger problem for people that include sequences that do not fully overlap.\n"; qanswers.push_back(qanswer);
        

        
        string issue = "Mothur can't find my input files. What wrong?"; issues.push_back(issue);
        string ianswer = "\tBy default, mothur will then look for the input files in the directory where mothur's executable is located. Mothur will also search the input, output and temporary default locations. You can set these locations using the set.dir command: set.dir(input=/users/myuser/desktop/mothurdata). Alternatively you can provide complete file names, or move the input files to mothur's executable location.\n"; ianswers.push_back(ianswer);
        
        issue = "I installed the latest version, but I am still running an older version. Why?"; issues.push_back(issue);
        ianswer = "\tWe often see this issue when you have an older version of mothur installed in your path. You can find out where by opening a terminal window and running: \n\n\tyourusername$ which mothur\n\tpath_to_old_version\n\tfor example: yourusername$ which mothur\n\t/usr/local/bin\n\n\tWhen you find the location of the older version, you can delete it or move it out of your path with the following:\n\n\tyourusername$ mv path_to_old_version/mothur new_location\n\tfor example: yourusername$ mv /usr/local/bin/mothur /Users/yourusername/desktop/old_version_mothur\n"; ianswers.push_back(ianswer);
        
        issue = "File Mismatches - 'yourSequence is in fileA but not in fileB, please correct.'"; issues.push_back(issue);
        ianswer = "\tThe most common reason this occurs is because you forgot to include a name or count file on a command, or accidentally included the wrong one due to a typo. Mothur has a 'current' option, which allows you to set file parameters to 'current'. For example, if fasta=current mothur will use the last fasta file given or created. The current option was designed to help avoid typo mistakes due to mothur's long filenames. Another reason this might occur is a process failing when you are using multiple processors. If a process dies, a file can be incomplete which would cause a mismatch error downstream.\n"; ianswers.push_back(ianswer);
        
        issue = "I don't have enough RAM or processing power. What are my options?"; issues.push_back(issue);
        ianswer = "\tIf you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required.\n\tAlternatively, you can use AWS to run your analysis. Here are instructions: https://mothur.org/wiki/Mothur_AMI.\n"; ianswers.push_back(ianswer);
        
        issue = "Mothur crashes when I read my distance file. What's wrong?"; issues.push_back(issue);
        ianswer = "\tThere are two common causes for this, file size and format.\n\n\tFileSize:\tThe cluster command loads your distance matrix into RAM, and your distance file is most likely too large to fit in RAM. There are two options to help with this. The first is to use a cutoff. By using a cutoff mothur will only load distances that are below the cutoff. If that is still not enough, there is a command called cluster.split, http://www.mothur.org/wiki/cluster.split. Cluster.split divides the dataset by taxonomic assignment and generates matrices for each grouping, and then clusters the smaller pieces separately. You may also be able to reduce the size of the original distance matrix by using the commands outline in the Schloss SOP, http://www.mothur.org/wiki/Schloss_SOP\n\n\tWrong Format:\tThis error can be caused by trying to read a column formatted distance matrix using the phylip parameter. By default, the dist.seqs command generates a column formatted distance matrix. To make a phylip formatted matrix set the dist.seqs command parameter output to lt.\n"; ianswers.push_back(ianswer);
        
        issue = "Why do I have such a large distance matrix?"; issues.push_back(issue);
        ianswer = "\tThis is most often caused by poor overlap of your reads. When reads have poor overlap, it greatly increases your error rate. Also, sequences that should cluster together don't because the errors appear to be genetic differences when in fact they are not. The quality of the data you are processing can not be overstressed. Error filled reads produce error filled results!\n\n\tCheck out Pat's blog: http://blog.mothur.org/2014/09/11/Why-such-a-large-distance-matrix/\n\n\tNOTE: To take a step back, if you look through our MiSeq SOP, you’ll see that we go to great pains to only work with the unique sequences to limit the number of sequences we have to align, screen for chimeras, classify, etc. We all know that 20 million reads will never make it through the pipeline without setting your computer on fire. Returning to the question at hand, you can imagine that if the reads do not fully overlap then any error in the 5’ end of the first read will be uncorrected by the 3’ end of the second read. If we assume for now that the errors are random, then every error will generate a new unique sequence. Granted, this happens less than 1% of the time, but multiply that by 20 million reads at whatever length you choose and you’ve got a big number. Viola, a bunch of unique reads and a ginormous distance matrix.\n"; ianswers.push_back(ianswer);
        
        issue = "Mothur reports a 'bad_alloc' error in the shhh.flows command. What's wrong?"; issues.push_back(issue);
        ianswer = "\tThis error indicates your computer is running out of memory. The shhh.flows command is very memory intensive. This error is most commonly caused by trying to process a dataset too large, using multiple processors, or failing to run trim.flows before shhh.flows. If you are using multiple processors, try running the command with processors=1, the more processors you use the more memory is required. Running trim.flows with an oligos file, and then shhh.flows with the file option may also resolve the issue. If for some reason you are unable to run shhh.flows with your data, a good alternative is to use the trim.seqs command using a 50-bp sliding window and to trim the sequence when the average quality score over that window drops below 35. Our results suggest that the sequencing error rates by this method are very good, but not quite as good as by shhh.flows and that the resulting sequences tend to be a bit shorter.\n"; ianswers.push_back(ianswer);
        
        
        string howto = "How do I make a tree?"; howtos.push_back(howto);
        string hanswer = "\tMothur has two commands that create trees: clearcut and tree.shared.\n\n\tThe clearcut commands creates a phylogenetic tree that represents how sequences relate. The clearcut program written by Initiative for Bioinformatics and Evolutionary Studies (IBEST) at the University of Idaho. For more information about clearcut please refer to http://bioinformatics.hungry.com/clearcut/\n\n\tThe tree.shared command will generate a newick-formatted tree file that describes the dissimilarity (1-similarity) among multiple groups. Groups are clustered using the UPGMA algorithm using the distance between communities as calculated using any of the calculators describing the similarity in community membership or structure.\n"; hanswers.push_back(hanswer);
        
        howto = "How do I know 'who' is in an OTU in a shared file?"; howtos.push_back(howto);
        hanswer = "\tYou can run the get.otulist command on the list file you used to generate the shared file. You want to be sure you are comparing the same distances. ie final.opti_mcc.0.03.otulist would relate to the 0.03 distance in your shared file. Also, if you subsample your data set and want to compare things, be sure to subsample the list and group file and then create the shared file to make sure you are working with the same sequences.\n\n\tsub.sample(list=yourListFile, count=yourCountFile, persample=t)\n\tmake.shared(list=yourSubsampledListFile, group=yourSubsampledCountFile, label=0.03)\n\tget.otulist(list=yourSubsampledListFile, label=0.03)\n"; hanswers.push_back(hanswer);
        
        howto = "How do I know 'who' is in the OTUs represented in the venn picture?"; howtos.push_back(howto);
        hanswer = "\tYou can use the get.sharedseqs command. Be sure to pay close attention to the 'unique' and 'shared' parameters.\n"; hanswers.push_back(hanswer);
        
        howto = "How do I select certain sequences or groups of sequences?"; howtos.push_back(howto);
        hanswer = "\tMothur has several 'get' and 'remove' commands: get.seqs, get.lineage, get.groups, get.dists, get.otus, remove.seqs, remove.lineage, remove.dists, remove.otus and remove.groups.\n"; hanswers.push_back(hanswer);
        
        howto = "How do I visualize my results from mothur?"; howtos.push_back(howto);
        hanswer = "\tTo visual your data with R follow this tutorial http://www.riffomonas.org/minimalR/06_line_plots.html.\n"; hanswers.push_back(hanswer);

        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "HelpCommand", "getCommonQuestions");
        exit(1);
    }
}
//**********************************************************************************************************************
int HelpCommand::execute(){
	try {
        if (commandName != "") {
            if (validCommands->isValidCommand(commandName)) {
                Command* command;
                string optionString = "help";
                
                if(commandName == "cluster")                    {	command = new ClusterCommand(optionString);					}
                else if(commandName == "unique.seqs")			{	command = new UniqueSeqsCommand(optionString);				}
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
                else if(commandName == "tree.shared")			{   command = new TreeSharedCommand(optionString);				}
                else if(commandName == "dist.shared")			{   command = new DistSharedCommand(optionString);			}
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
                else if((commandName == "get.otus")	|| (commandName == "get.otulabels"))			{	command = new GetOtusCommand(optionString);			}
                else if((commandName == "remove.otus") || (commandName == "remove.otulabels"))			{	command = new RemoveOtusCommand(optionString);			}
                else if((commandName == "list.otus")	||(commandName == "list.otulabels"))        {	command = new ListOtusCommand(optionString);           }
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
                m->mothurOut("[ERROR]: " + commandName + " is not a valid command.\n");
                validCommands->printCommands(cout);
            }
        }else {
            //validCommands->printCommands(cout);
            validCommands->printCommandsCategories(cout);
#if defined NON_WINDOWS
            cout << BOLDMAGENTA << "\nFor more information about a specific command type 'commandName(help)' i.e. 'cluster(help)'\n"; cout << RESET << endl;
            m->mothurOutJustToLog("\nFor more information about a specific command type 'commandName(help)' i.e. 'cluster(help)'\n");
#else
            m->mothurOut("\nFor more information about a specific command type 'commandName(help)' i.e. 'cluster(help)'\n");
#endif
            getCommonQuestions();
        }
		
#if defined NON_WINDOWS
        cout << BOLDMAGENTA << "\nFor further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki, or contact Pat Schloss at mothur.bugs@gmail.com.\n"; cout << RESET << endl;
        m->mothurOutJustToLog("\nFor further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki, or contact Pat Schloss at mothur.bugs@gmail.com.\n");
#else
        m->mothurOut("\nFor further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki, or contact Pat Schloss at mothur.bugs@gmail.com.\n");
#endif
        
	
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HelpCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************/
