#include "distcdataset.h"
#include "getdistscommand.h"
#include "listseqscommand.h"
#include "getseqscommand.h"

/***********************************************************************/
DistCDataSet::DistCDataSet() {
    m = MothurOut::getInstance();
    columnFile = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/stability.MISeq_SOP.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.dist";
    countFile = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/stability.count_table";
}
/***********************************************************************/
vector<string> DistCDataSet::getFiles(int numSeqs) {
    vector<string> newFiles;
    
    if (numSeqs > 2055) { m->mothurOut("[ERROR]: too many seqs requested in DistCDataSet::getFiles\n"); }
    else {
        string inputString = "count=" + countFile;
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: list.seqs(" + inputString + ")"); m->mothurOutEndLine();
        m->setMothurCalling(true);
        
        Command* listCommand = new ListSeqsCommand(inputString);
        listCommand->execute();
        
        map<string, vector<string> > filenames = listCommand->getOutputFiles();
        
        delete listCommand;
        m->setMothurCalling(false);
        
        string accnosfile = filenames["accnos"][0];
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        ifstream in;
        m->openInputFile(accnosfile, in);
        
        ofstream out;
        m->openOutputFile("temp.accnos", out);
        
        int count = 0; string name;
        while(!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            in >> name; m->gobble(in);
            out << name << endl;
            count++;
            
            if (count >= numSeqs) { break; }
        }
        in.close();
        out.close();
        m->mothurRemove(accnosfile);
        
        inputString = "count=" + countFile + ", accnos=temp.accnos";
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: get.seqs(" + inputString + ")"); m->mothurOutEndLine();
        m->setMothurCalling(true);
        
        Command* getCommand = new GetSeqsCommand(inputString);
        getCommand->execute();
        
        filenames = getCommand->getOutputFiles();
        
        delete getCommand;
        m->setMothurCalling(false);
        
        string newCountfile = filenames["count"][0];
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        inputString = "column=" + columnFile + ", accnos=temp.accnos";
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: get.dists(" + inputString + ")"); m->mothurOutEndLine();
        m->setMothurCalling(true);
        
        Command* getDCommand = new GetDistsCommand(inputString);
        getDCommand->execute();
        
        filenames = getDCommand->getOutputFiles();
        
        delete getDCommand;
        m->setMothurCalling(false);
        
        string newColumnfile = filenames["column"][0];
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        newFiles.push_back(newColumnfile); newFiles.push_back(newCountfile);
    }
    
    return newFiles;
}
/***********************************************************************/
