/*
 *  aligncommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 *	This version of nast does everything I think that the greengenes nast server does and then some.  I have added the 
 *	feature of allowing users to define their database, kmer size for searching, alignment penalty values and alignment 
 *	method.  This latter feature is perhaps most significant.  nastPlus enables a user to use either a Needleman-Wunsch 
 *	(non-affine gap penalty) or Gotoh (affine gap penalty) pairwise alignment algorithm.  This is significant because it
 *	allows for a global alignment and not the local alignment provided by bLAst.  Furthermore, it has the potential to
 *	provide a better alignment because of the banding method employed by blast (I'm not sure about this).
 *
 */

#include "aligncommand.h"

//**********************************************************************************************************************
vector<string> AlignCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pcandidate("fasta", "InputTypes", "", "", "none", "none", "none","fasta-alignreport-accnos",false,true,true); parameters.push_back(pcandidate);
		CommandParameter psearch("search", "Multiple", "kmer-blast-suffix", "kmer", "", "", "","",false,false,true); parameters.push_back(psearch);
		CommandParameter pksize("ksize", "Number", "", "8", "", "", "","",false,false); parameters.push_back(pksize);
		CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
		CommandParameter palign("align", "Multiple", "needleman-gotoh-blast-noalign", "needleman", "", "", "","",false,false,true); parameters.push_back(palign);
		CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
		CommandParameter pgapopen("gapopen", "Number", "", "-5.0", "", "", "","",false,false); parameters.push_back(pgapopen);
		CommandParameter pgapextend("gapextend", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapextend);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pflip("flip", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pflip);
		CommandParameter pthreshold("threshold", "Number", "", "0.50", "", "", "","",false,false); parameters.push_back(pthreshold);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["alignreport"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;

        abort = false; calledHelp = false;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string AlignCommand::getHelpString(){	
	try {
		string helpString = "\n";
		helpString += "The align.seqs command reads a file containing sequences and creates an alignment file and a report file.\n";
		helpString += "The align.seqs command parameters are " + getCommandParameters() + ".\n";
		helpString += "The reference and fasta parameters are required. You may leave fasta blank if you have a valid fasta file.\n";
		helpString += "The search parameter allows you to specify the method to find most similar reference sequence.  Your options are: suffix, kmer and blast. The default is kmer.\n";
		helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n";
		helpString += "The ksize parameter allows you to specify the kmer size for finding most similar reference to a given sequence.  The default is 8.\n";
		helpString += "The match parameter allows you to specify the bonus for having the same base. Default=1.0.\n";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases. Default=-1.0.\n";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. Default=-5.0.\n";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment. Default=-2.0.\n";
        helpString += "If the flip parameter is set to true the reverse complement of the sequence is aligned and the better alignment is reported.";
		helpString += " By default, mothur will align the reverse compliment of your sequences when the alignment process removes more than 50% of the bases indicating the read may be flipped. This process assembles the best possible alignment, and downstream analysis will remove any poor quality reads remaining.\n";
		helpString += "The threshold is used to specify a cutoff at which an alignment is deemed 'bad' and the reverse complement may be tried. The default threshold is 0.50, meaning 50% of the bases are removed in the alignment.\n";
		helpString += "The align.seqs command should be in the following format: ";
		helpString += "align.seqs(reference=yourTemplateFile, fasta=yourUnalignedFastaFile)\n";
		helpString += "Example: align.seqs(fasta=water.fasta, template=silva.v4.fasta)\n\n";
        
        getCommonQuestions();

		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string AlignCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string issue = "...template is not aligned, aborting. What do I do?"; issues.push_back(issue);
        string ianswer = "\tMothur requires the reference file to be aligned to generate aligned sequences. You can download mothur's aligned silva references here, https://mothur.org/wiki/Silva_reference_files. For ITS sequences, see 'how to' below.\n"; ianswers.push_back(ianswer);
        
        issue = "...xxx of your sequences generated alignments that eliminated too many bases... What does this mean?"; issues.push_back(issue);
        ianswer = "\tBy default, mothur will align the reverse compliment of your sequences when the alignment process removes more than 50% of the bases indicating the read may be flipped. This process assembles the best possible alignment, and downstream analysis will remove any poor quality reads remaining.\n"; ianswers.push_back(ianswer);
        
        
        string howto = "How do I 'align' ITS sequences?"; howtos.push_back(howto);
        string hanswer = "\tYou really can't do an alignment because there isn't positional homology. You can use the pre.cluster and pairwise.seqs commands to generate a distance matrix from unaligned sequences.\n"; hanswers.push_back(hanswer);
        
        howto = "How do I create a custom reference for the region I am studying?"; howtos.push_back(howto);
        hanswer = "\tYou can tailor your reference using this method: http://blog.mothur.org/2016/07/07/Customization-for-your-region/.\n"; hanswers.push_back(hanswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);

        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "AlignCommand", "getCommonQuestions");
        exit(1);
    }
}
//**********************************************************************************************************************
string AlignCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],align"; } //makes file like: amazon.align
        else if (type == "alignreport") {  pattern = "[filename],align.report"; }
        else if (type == "accnos") {  pattern = "[filename],flip.accnos"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "AlignCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
AlignCommand::AlignCommand(string option)  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true;}
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
            ValidParameters validParameter;
            templateFileName = validParameter.validFile(parameters, "reference");
            if (templateFileName == "not found") { m->mothurOut("[ERROR]: The reference parameter is a required for the align.seqs command, aborting.\n"); abort = true;
            }else if (templateFileName == "not open") { abort = true; }
            
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not found") {
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
                else { 	m->mothurOut("[ERROR]: You have no current fasta file and the fasta parameter is required.\n");  abort = true; }
            }
            else if (fastafile == "not open") { abort = true; }
            else { current->setFastaFile(fastafile); }
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.valid(parameters, "ksize");		if (temp == "not found"){	temp = "8";				}
			util.mothurConvert(temp, kmerSize); 
			
			temp = validParameter.valid(parameters, "match");		if (temp == "not found"){	temp = "1.0";			}
			util.mothurConvert(temp, match);  
			
			temp = validParameter.valid(parameters, "mismatch");		if (temp == "not found"){	temp = "-1.0";			}
			util.mothurConvert(temp, misMatch);  
			
			temp = validParameter.valid(parameters, "gapopen");		if (temp == "not found"){	temp = "-5.0";			}
			util.mothurConvert(temp, gapOpen);  
			
			temp = validParameter.valid(parameters, "gapextend");	if (temp == "not found"){	temp = "-2.0";			}
			util.mothurConvert(temp, gapExtend); 
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
			temp = validParameter.valid(parameters, "flip");			if (temp == "not found"){	temp = "t";				}
			flip = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "threshold");	if (temp == "not found"){	temp = "0.50";			}
			util.mothurConvert(temp, threshold); 
			
			search = validParameter.valid(parameters, "search");		if (search == "not found"){	search = "kmer";		}
			if ((search != "suffix") && (search != "kmer") && (search != "blast")) { m->mothurOut("invalid search option: choices are kmer, suffix or blast.\n");  abort=true; }
			
			align = validParameter.valid(parameters, "align");		if (align == "not found"){	align = "needleman";	}
			if ((align != "needleman") && (align != "gotoh") && (align != "blast") && (align != "noalign")) { m->mothurOut("invalid align option: choices are needleman, gotoh, blast or noalign.\n");  abort=true; }
            
            if ((search == "blast") || (align == "blast")) {
                string blastlocation = "";
                vector<string> locations = current->getLocations();
                string path = current->getProgramPath();
                
                if (!util.findBlastLocation(blastlocation, path, locations)){
                    m->mothurOut("[WARNING]: Unable to locate blast executables, cannot use blast as search or align method. Using defaults instead.\n");
                    if (search == "blast") { search = "kmer"; }
                    if (align == "blast") { align = "needleman"; }
                }
            }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "AlignCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
AlignCommand::~AlignCommand(){}
//**********************************************************************************************************************

int AlignCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        long long before = util.getRAMUsed(); long long total = util.getTotalRAM();
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: RAM used before reading template " + toString(before) + " of total RAM available " + toString(total) + "\n"); }
        
        templateDB = new AlignmentDB(templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch, util.getRandomNumber(), true);
        
        if (m->getControl_pressed()) { outputTypes.clear(); return 0; }
        
        time_t start = time(NULL);
        m->mothurOut("\nAligning sequences from " + fastafile + " ...\n" );
        
        if (outputdir == "") {  outputdir += util.hasPath(fastafile); }
        map<string, string> variables; variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        string alignFileName = getOutputFileName("fasta", variables);
        string reportFileName = getOutputFileName("alignreport", variables);
        string accnosFileName = getOutputFileName("accnos", variables);
        
        bool hasAccnos = true;
        vector<long long> numFlipped;
        numFlipped.push_back(0); //numflipped because reverse was better
        numFlipped.push_back(0); //total number of sequences with over 50% of bases removed
        
        long long numFastaSeqs = createProcesses(alignFileName, reportFileName, accnosFileName, fastafile, numFlipped);
        
        delete templateDB;
        
        if (m->getControl_pressed()) { util.mothurRemove(accnosFileName); util.mothurRemove(alignFileName); util.mothurRemove(reportFileName); outputTypes.clear();  return 0; }
        
        //delete accnos file if its blank else report to user
        if (util.isBlank(accnosFileName)) {  util.mothurRemove(accnosFileName);  hasAccnos = false; }
        else {
            m->mothurOut("[WARNING]: " + toString(numFlipped[1]) + " of your sequences generated alignments that eliminated too many bases, a list is provided in " + accnosFileName + ".");
            if (!flip) {
                m->mothurOut(" If you set the flip parameter to true mothur will try aligning the reverse compliment as well. flip=t");
            }else{  m->mothurOut("\n[NOTE]: " + toString(numFlipped[0]) + " of your sequences were reversed to produce a better alignment.");  }
            m->mothurOutEndLine();
        }
        
        outputNames.push_back(alignFileName); outputTypes["fasta"].push_back(alignFileName);
        outputNames.push_back(reportFileName); outputTypes["alignreport"].push_back(reportFileName);
        if (hasAccnos)	{	outputNames.push_back(accnosFileName);	outputTypes["accnos"].push_back(accnosFileName);  }
		
        m->mothurOut("\nIt took " + toString(time(NULL) - start) + " seconds to align " + toString(numFastaSeqs) + " sequences.\n");
        
		//set align file as new current fastafile
		string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; current->setFastaFile(currentFasta); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
struct alignStruct {
    OutputWriter* alignWriter;
    OutputWriter* reportWriter;
    OutputWriter* accnosWriter;
    string inputFilename;
    string alignMethod, search, templateFileName;
    float match, misMatch, gapOpen, gapExtend, threshold;
    bool flip;
    long long numSeqs;
    int kmerSize;
    
    vector<long long> flippedResults;
    linePair filePos;
    
    MothurOut* m;
    Utils util;
    AlignmentDB* templateDB;
    Alignment* alignment;
    
    alignStruct (linePair fP, OutputWriter* aFName, OutputWriter* reFName, OutputWriter* ac, string fname, AlignmentDB* tfn, string al, float ma, float misMa, float gOpen, float gExtend, float thr, bool fl, int ks, string se) {
        
        filePos.start = fP.start;
        filePos.end = fP.end;
        alignWriter = aFName;
        reportWriter = reFName;
        accnosWriter = ac;
        inputFilename = fname;
        numSeqs = 0;
        m = MothurOut::getInstance();
        match = ma;
        misMatch = misMa;
        gapOpen = gOpen;
        gapExtend = gExtend;
        threshold = thr;
        flip = fl;
        search = se;
        kmerSize = ks;
        flippedResults.resize(2, 0);
        
        templateDB = tfn;
        
        int longestBase = templateDB->getLongestBase();
        if (m->getDebug()) { m->mothurOut("[DEBUG]: template longest base = "  + toString(longestBase) + " \n");            }
        if(al == "gotoh")            {    alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);   }
        else if(al == "needleman")    {    alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);         }
        else if(al == "blast")        {    alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);             }
        else if(al == "noalign")        {    alignment = new NoAlign();                                                     }
        else {
            m->mothurOut(al + " is not a valid alignment option. I will run the command using needleman.\n");
            alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);
        }
    }
    ~alignStruct() { delete alignment;  }
    
};
//**********************************************************************************************************************
void alignDriver(alignStruct* params) {
	try {
        NastReport report;
		
		ifstream inFASTA;
		params->util.openInputFile(params->inputFilename, inFASTA);

		inFASTA.seekg(params->filePos.start);

		bool done = false;
        
		long long count = 0;
        long long numFlipped_0 = 0;
        long long numFlipped_1 = 0;
		
		while (!done) {
			
			if (params->m->getControl_pressed()) {  break; }
			
			Sequence* candidateSeq = new Sequence(inFASTA);  params->util.gobble(inFASTA);
			report.setCandidate(candidateSeq);

			int origNumBases = candidateSeq->getNumBases();
			string originalUnaligned = candidateSeq->getUnaligned();
			int numBasesNeeded = origNumBases * params->threshold;
	
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getUnaligned().length()+1 > params->alignment->getnRows()) {
                    if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + candidateSeq->getName() + " " + toString(candidateSeq->getUnaligned().length()) + " " + toString(params->alignment->getnRows()) + " \n"); }
					params->alignment->resize(candidateSeq->getUnaligned().length()+2);
				}
                
                float searchScore;
				Sequence temp = params->templateDB->findClosestSequence(candidateSeq, searchScore);
				Sequence* templateSeq = new Sequence(temp.getName(), temp.getAligned());
								
				Nast* nast = new Nast(params->alignment, candidateSeq, templateSeq);
		
				Sequence* copy;
				
				Nast* nast2;
				bool needToDeleteCopy = false;  //this is needed in case you have you enter the ifs below
												//since nast does not make a copy of hte sequence passed, and it is used by the reporter below
												//you can't delete the copy sequence til after you report, but you may choose not to create it in the first place
												//so this bool tells you if you need to delete it
												
				//if there is a possibility that this sequence should be reversed
				if (candidateSeq->getNumBases() < numBasesNeeded) {
					numFlipped_1++;
					string wasBetter =  "";
					//if the user wants you to try the reverse
					if (params->flip) {
				
						//get reverse compliment
						copy = new Sequence(candidateSeq->getName(), originalUnaligned);
						copy->reverseComplement();
                        
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: flipping "  + candidateSeq->getName() + " \n"); }
						
						//rerun alignment
						Sequence temp2 = params->templateDB->findClosestSequence(copy, searchScore);
						Sequence* templateSeq2 = new Sequence(temp2.getName(), temp2.getAligned());
                        
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: closest template "  + temp2.getName() + " \n"); }
						
						nast2 = new Nast(params->alignment, copy, templateSeq2);
                        
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: completed Nast2 "  + candidateSeq->getName() + " flipped numBases = " + toString(copy->getNumBases()) + " old numbases = " + toString(candidateSeq->getNumBases()) +" \n"); }
			
						//check if any better
						if (copy->getNumBases() > candidateSeq->getNumBases()) {
							candidateSeq->setAligned(copy->getAligned());  //use reverse compliments alignment since its better
                            delete templateSeq;
							templateSeq = templateSeq2;
							delete nast;
							nast = nast2;
							needToDeleteCopy = true;
							wasBetter = "\treverse complement produced a better alignment, so mothur used the reverse complement.";
                            numFlipped_0++;
						}else{  
							wasBetter = "\treverse complement did NOT produce a better alignment so it was not used, please check sequence.";
							delete nast2;
                            delete templateSeq2;
							delete copy;	
						}
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: done.\n"); }
					}
					
					//create accnos file with names
					params->accnosWriter->write(candidateSeq->getName() + wasBetter + "\n");
				}
				
				report.setTemplate(templateSeq);
				report.setSearchParameters(params->search, searchScore);
				report.setAlignmentParameters(params->alignMethod, params->alignment);
				report.setNastParameters(*nast);
	
				params->alignWriter->write('>' + candidateSeq->getName() + '\n' + candidateSeq->getAligned() + '\n');
				params->reportWriter->write(report.getReport());
				delete nast;
                delete templateSeq;
				if (needToDeleteCopy) {   delete copy;   }
                
				count++;
			}
			delete candidateSeq;
			
			#if defined NON_WINDOWS
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= params->filePos.end)) { break; }
			#else
				if (count == params->filePos.end) { break; }
			#endif
			
			//report progress
			if((count) % 1000 == 0){	params->m->mothurOutJustToScreen(toString(count) + "\n"); 		}
			
		}
		//report progress
		if((count) % 1000 != 0){	params->m->mothurOutJustToScreen(toString(count) + "\n"); 		}
        
        params->numSeqs += count;
        params->flippedResults[0] += numFlipped_0;
        params->flippedResults[1] += numFlipped_1;
        
		inFASTA.close();
		
	}
	catch(exception& e) {
		params->m->errorOut(e, "AlignCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
//void alignDriver(linePair* filePos, string alignFName, string reportFName, string accnosFName, string filename, vector<long long>& numFlipped,MothurOut* m, string align, float match, float misMatch, float gapOpen, float gapExtend, float threshold, bool flip, AlignmentDB* templateDB, string search, long long& count) {
long long AlignCommand::createProcesses(string alignFileName, string reportFileName, string accnosFName, string filename, vector<long long>& numFlipped) {
	try {
        
        vector<linePair> lines;
        vector<double> positions;
#if defined NON_WINDOWS
        positions = util.divideFile(filename, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        long long numFastaSeqs = 0;
        positions = util.setFilePosFasta(filename, numFastaSeqs);
        if (numFastaSeqs < processors) { processors = numFastaSeqs; m->mothurOut("Reducing processors to " + toString(numFastaSeqs) + ".\n");  }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numFastaSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<alignStruct*> data;
        
        long long num = 0;
        for (int i = 0; i < numFlipped.size(); i++) { numFlipped[i] = 0; }
        
        time_t start, end;
        time(&start);
        
        NastReport nast; string nastHeaders = nast.getHeaders();
        ofstream out; util.openOutputFile(reportFileName, out); out << nastHeaders; out.close();
        
        auto synchronizedOutputAlignFile = std::make_shared<SynchronizedOutputFile>(alignFileName);
        auto synchronizedOutputReportFile = std::make_shared<SynchronizedOutputFile>(reportFileName, true);
        auto synchronizedOutputAccnosFile = std::make_shared<SynchronizedOutputFile>(accnosFName);

        for (int i = 0; i < processors-1; i++) {
            
            OutputWriter* threadAlignWriter = new OutputWriter(synchronizedOutputAlignFile);
            OutputWriter* threadReportWriter = new OutputWriter(synchronizedOutputReportFile);
            OutputWriter* threadAccnosWriter = new OutputWriter(synchronizedOutputAccnosFile);

            
            alignStruct* dataBundle = new alignStruct(lines[i+1], threadAlignWriter, threadReportWriter, threadAccnosWriter, filename, templateDB, align, match, misMatch, gapOpen, gapExtend, threshold, flip, kmerSize, search);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(alignDriver, dataBundle));
         }
        
        OutputWriter* threadAlignWriter = new OutputWriter(synchronizedOutputAlignFile);
        OutputWriter* threadReportWriter = new OutputWriter(synchronizedOutputReportFile);
        OutputWriter* threadAccnosWriter = new OutputWriter(synchronizedOutputAccnosFile);
        
        alignStruct* dataBundle = new alignStruct(lines[0], threadAlignWriter, threadReportWriter, threadAccnosWriter, filename, templateDB, align, match, misMatch, gapOpen, gapExtend, threshold, flip, kmerSize, search);
        alignDriver(dataBundle);
        numFlipped[0] = dataBundle->flippedResults[0];
        numFlipped[1] = dataBundle->flippedResults[1];
        num = dataBundle->numSeqs;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->numSeqs;
            numFlipped[0] += data[i]->flippedResults[0];
            numFlipped[1] += data[i]->flippedResults[1];
            
            delete data[i]->alignWriter;
            delete data[i]->reportWriter;
            delete data[i]->accnosWriter;
            
            delete data[i];
            delete workerThreads[i];
        }
        synchronizedOutputAlignFile->close(); synchronizedOutputReportFile->close(); synchronizedOutputAccnosFile->close();
        delete threadAlignWriter; delete threadAccnosWriter; delete threadReportWriter;
        delete dataBundle;
        
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to align " + toString(num) + " sequences.\n\n");
        
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************
