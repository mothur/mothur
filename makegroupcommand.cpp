/*
 *  makegroupcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/7/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "makegroupcommand.h"

//**********************************************************************************************************************

MakeGroupCommand::MakeGroupCommand(string option)  {
	try {
		
		abort = false;
	
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"fasta","groups","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}

			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			

			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			

			fastaFileName = validParameter.validFile(parameters, "fasta", false);
			if (fastaFileName == "not found") { m->mothurOut("fasta is a required parameter for the make.group command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					if (inputDir != "") {
						string path = hasPath(fastaFileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
					}
	
					int ableToOpen;
					ifstream in;
					
					ableToOpen = openInputFile(fastaFileNames[i], in);
					in.close();

					if (ableToOpen == 1) { 
						m->mothurOut(fastaFileNames[i] + " will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						fastaFileNames.erase(fastaFileNames.begin()+i);
						i--;
					}
					
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
		
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "MakeGroupCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

MakeGroupCommand::~MakeGroupCommand(){	}

//**********************************************************************************************************************

void MakeGroupCommand::help(){
	try {
		m->mothurOut("The align.seqs command reads a file containing sequences and creates an alignment file and a report file.\n");
		m->mothurOut("The align.seqs command parameters are template, candidate, search, ksize, align, match, mismatch, gapopen and gapextend.\n");
		m->mothurOut("The template and candidate parameters are required. You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n");
		m->mothurOut("The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer.\n");
		m->mothurOut("The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n");
		m->mothurOut("The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 8.\n");
		m->mothurOut("The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n");
		m->mothurOut("The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n");
		m->mothurOut("The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n");
		m->mothurOut("The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n");
		m->mothurOut("The flip parameter is used to specify whether or not you want mothur to try the reverse complement if a sequence falls below the threshold.  The default is false.\n");
		m->mothurOut("The threshold is used to specify a cutoff at which an alignment is deemed 'bad' and the reverse complement may be tried. The default threshold is 0.50, meaning 50% of the bases are removed in the alignment.\n");
		m->mothurOut("If the flip parameter is set to true the reverse complement of the sequence is aligned and the better alignment is reported.\n");
		m->mothurOut("The default for the threshold parameter is 0.50, meaning at least 50% of the bases must remain or the sequence is reported as potentially reversed.\n");
		m->mothurOut("The align.seqs command should be in the following format: \n");
		m->mothurOut("align.seqs(template=yourTemplateFile, candidate=yourCandidateFile, align=yourAlignmentMethod, search=yourSearchmethod, ksize=yourKmerSize, match=yourMatchBonus, mismatch=yourMismatchpenalty, gapopen=yourGapopenPenalty, gapextend=yourGapExtendPenalty) \n");
		m->mothurOut("Example align.seqs(candidate=candidate.fasta, template=core.filtered, align=kmer, search=gotoh, ksize=8, match=2.0, mismatch=3.0, gapopen=-2.0, gapextend=-1.0)\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. candidate), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

int MakeGroupCommand::execute(){
	try {
		if (abort == true) {	return 0;	}
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		//for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "execute");
		exit(1);
	}
}

