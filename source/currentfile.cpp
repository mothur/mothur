//
//  currentfile.cpp
//  Mothur
//
//  Created by Sarah Westcott on 11/9/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "currentfile.h"

/*********************************************************************************************/
set<string> CurrentFile::getCurrentTypes()  {
    try {
        
        set<string> types;
        types.insert("fasta");
        types.insert("summary");
        types.insert("file");
        types.insert("accnos");
        types.insert("column");
        types.insert("design");
        types.insert("group");
        types.insert("list");
        types.insert("name");
        types.insert("oligos");
        types.insert("order");
        types.insert("ordergroup");
        types.insert("phylip");
        types.insert("qfile");
        types.insert("relabund");
        types.insert("clr");
        types.insert("sabund");
        types.insert("rabund");
        types.insert("sff");
        types.insert("shared");
        types.insert("taxonomy");
        types.insert("constaxonomy");
        types.insert("contigsreport");
        types.insert("tree");
        types.insert("flow");
        types.insert("biom");
        types.insert("count");
        types.insert("processors");
        types.insert("sample");
        
        return types;
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "getCurrentTypes");
        exit(1);
    }
}
/*********************************************************************************************/
void CurrentFile::printCurrentFiles(string filename)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        if (filename != "") {
            ofstream out; util.openOutputFile(filename, out);
            
            if (accnosfile != "")		{  out << "accnos=" + accnosfile + "\n";           }
            if (columnfile != "")		{  out << "column=" + columnfile + "\n";			}
            if (designfile != "")		{  out << "design=" + designfile + "\n";		}
            if (fastafile != "")		{  out << "fasta=" + fastafile + "\n";			}
            if (groupfile != "")		{  out << "group=" + groupfile + "\n";			}
            if (listfile != "")			{  out << "list=" + listfile + "\n";				}
            if (namefile != "")			{  out << "name=" + namefile + "\n";			}
            if (oligosfile != "")		{  out << "oligos=" + oligosfile + "\n";			}
            if (orderfile != "")		{  out << "order=" + orderfile + "\n";			}
            if (ordergroupfile != "")	{  out << "ordergroup=" + ordergroupfile + "\n";	}
            if (phylipfile != "")		{  out << "phylip=" + phylipfile + "\n";		}
            if (qualfile != "")			{  out << "qfile=" + qualfile + "\n";				}
            if (rabundfile != "")		{  out << "rabund=" + rabundfile + "\n";			}
            if (relabundfile != "")		{  out << "relabund=" + relabundfile + "\n";		}
            if (clrfile != "")          {  out << "clr=" + clrfile + "\n";               }
            if (sabundfile != "")		{  out << "sabund=" + sabundfile + "\n";			}
            if (sfffile != "")			{  out << "sff=" + sfffile + "\n";			}
            if (sharedfile != "")		{  out << "shared=" + sharedfile + "\n";			}
            if (taxonomyfile != "")		{  out << "taxonomy=" + taxonomyfile + "\n";		}
            if (constaxonomyfile != "")	{  out << "constaxonomy=" + constaxonomyfile + "\n"; }
            if (contigsreportfile != ""){  out << "contigsreport=" + contigsreportfile + "\n";}
            if (treefile != "")			{  out << "tree=" + treefile + "\n";			}
            if (flowfile != "")			{  out << "flow=" + flowfile + "\n";			}
            if (biomfile != "")			{  out << "biom=" + biomfile + "\n";			}
            if (countfile != "")        {  out << "count=" + countfile + "\n";      }
            if (processors != "1")		{  out << "processors=" + processors + "\n";		}
            if (summaryfile != "")		{  out << "summary=" + summaryfile + "\n";       }
            if (filefile != "")         {  out << "file=" + filefile + "\n";           }
            if (samplefile != "")       {  out << "sample=" + samplefile + "\n";              }
            
            out.close();
        }
        
        if (accnosfile != "")		{  m->mothurOut("accnos=" + accnosfile); m->mothurOutEndLine();			}
        if (columnfile != "")		{  m->mothurOut("column=" + columnfile); m->mothurOutEndLine();			}
        if (designfile != "")		{  m->mothurOut("design=" + designfile); m->mothurOutEndLine();			}
        if (fastafile != "")		{  m->mothurOut("fasta=" + fastafile); m->mothurOutEndLine();				}
        if (groupfile != "")		{  m->mothurOut("group=" + groupfile); m->mothurOutEndLine();				}
        if (listfile != "")			{  m->mothurOut("list=" + listfile); m->mothurOutEndLine();				}
        if (namefile != "")			{  m->mothurOut("name=" + namefile); m->mothurOutEndLine();				}
        if (oligosfile != "")		{  m->mothurOut("oligos=" + oligosfile); m->mothurOutEndLine();			}
        if (orderfile != "")		{  m->mothurOut("order=" + orderfile); m->mothurOutEndLine();				}
        if (ordergroupfile != "")	{  m->mothurOut("ordergroup=" + ordergroupfile); m->mothurOutEndLine();	}
        if (phylipfile != "")		{  m->mothurOut("phylip=" + phylipfile); m->mothurOutEndLine();			}
        if (qualfile != "")			{  m->mothurOut("qfile=" + qualfile); m->mothurOutEndLine();				}
        if (rabundfile != "")		{  m->mothurOut("rabund=" + rabundfile); m->mothurOutEndLine();			}
        if (relabundfile != "")		{  m->mothurOut("relabund=" + relabundfile); m->mothurOutEndLine();		}
        if (clrfile != "")          {  m->mothurOut("clr=" + clrfile); m->mothurOutEndLine();               }
        if (sabundfile != "")		{  m->mothurOut("sabund=" + sabundfile); m->mothurOutEndLine();			}
        if (sfffile != "")			{  m->mothurOut("sff=" + sfffile); m->mothurOutEndLine();					}
        if (sharedfile != "")		{  m->mothurOut("shared=" + sharedfile); m->mothurOutEndLine();			}
        if (taxonomyfile != "")		{  m->mothurOut("taxonomy=" + taxonomyfile); m->mothurOutEndLine();		}
        if (constaxonomyfile != "")	{  m->mothurOut("constaxonomy=" + constaxonomyfile); m->mothurOutEndLine();}
        if (contigsreportfile != ""){  m->mothurOut("contigsreport=" + contigsreportfile); m->mothurOutEndLine();}
        if (treefile != "")			{  m->mothurOut("tree=" + treefile); m->mothurOutEndLine();				}
        if (flowfile != "")			{  m->mothurOut("flow=" + flowfile); m->mothurOutEndLine();				}
        if (biomfile != "")			{  m->mothurOut("biom=" + biomfile); m->mothurOutEndLine();				}
        if (countfile != "")        {  m->mothurOut("count=" + countfile); m->mothurOutEndLine();        }
        if (processors != "1")		{  m->mothurOut("processors=" + processors); m->mothurOutEndLine();		}
        if (summaryfile != "")		{  m->mothurOut("summary=" + summaryfile); m->mothurOutEndLine();         }
        if (filefile != "")         {  m->mothurOut("file=" + filefile); m->mothurOutEndLine();               }
        if (samplefile != "")       {  m->mothurOut("sample=" + samplefile); m->mothurOutEndLine();               }
        
        
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "printCurrentFiles");
        exit(1);
    }
}
/*********************************************************************************************/
bool CurrentFile::hasCurrentFiles()  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        bool hasCurrent = false;
        
        if (accnosfile != "")		{  return true;			}
        if (columnfile != "")		{  return true;			}
        if (designfile != "")		{  return true;			}
        if (fastafile != "")		{  return true;			}
        if (groupfile != "")		{  return true;			}
        if (listfile != "")			{  return true;			}
        if (namefile != "")			{  return true;			}
        if (oligosfile != "")		{  return true;			}
        if (orderfile != "")		{  return true;			}
        if (ordergroupfile != "")	{  return true;			}
        if (phylipfile != "")		{  return true;			}
        if (qualfile != "")			{  return true;			}
        if (rabundfile != "")		{  return true;			}
        if (relabundfile != "")		{  return true;			}
        if (clrfile != "")          {  return true;         }
        if (sabundfile != "")		{  return true;			}
        if (sfffile != "")			{  return true;			}
        if (sharedfile != "")		{  return true;			}
        if (taxonomyfile != "")		{  return true;			}
        if (constaxonomyfile != "")	{  return true;			}
        if (contigsreportfile != ""){  return true;			}
        if (treefile != "")			{  return true;			}
        if (flowfile != "")			{  return true;			}
        if (biomfile != "")			{  return true;			}
        if (countfile != "")        {  return true;			}
        if (summaryfile != "")      {  return true;			}
        if (filefile != "")         {  return true;			}
        if (samplefile != "")       {  return true;         }
        if (processors != "1")		{  return true;			}
        
        return hasCurrent;
        
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "hasCurrentFiles");
        exit(1);
    }
}

/*********************************************************************************************/
void CurrentFile::clearCurrentFiles()  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        phylipfile = "";
        filefile = "";
        columnfile = "";
        listfile = "";
        rabundfile = "";
        sabundfile = "";
        namefile = "";
        groupfile = "";
        designfile = "";
        orderfile = "";
        treefile = "";
        sharedfile = "";
        ordergroupfile = "";
        contigsreportfile = "";
        constaxonomyfile = "";
        relabundfile = "";
        clrfile = "";
        fastafile = "";
        qualfile = "";
        sfffile = "";
        oligosfile = "";
        accnosfile = "";
        taxonomyfile = "";
        flowfile = "";
        biomfile = "";
        countfile = "";
        summaryfile = "";
        samplefile = "";
        unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
        if (concurentThreadsSupported < 1) { concurentThreadsSupported = 1; } //in case thread errors
        processors = toString(concurentThreadsSupported);
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "clearCurrentFiles");
        exit(1);
    }
}
/*********************************************************************************************/
int CurrentFile::setProcessors(string p)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        if (!util.isInteger(p)) {
            unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
            if (concurentThreadsSupported < 1) { concurentThreadsSupported = 1; } //in case thread errors
            processors = toString(concurentThreadsSupported);
            m->mothurOut("[ERROR]: " + p + " is not an integer. Setting processors to " + toString(processors) + "\n");
        }else { processors = p;  m->mothurOut("\nUsing " + toString(processors) + " processors.\n"); }
        int numProcessors = 1;
        util.mothurConvert(p, numProcessors);
        return numProcessors;
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "clearCurrentFiles");
        exit(1);
    }
}
/*********************************************************************************************/
void CurrentFile::setDefaultPath(string pathname)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        if (pathname != "") { //add / to name if needed
            string lastChar = pathname.substr(pathname.length()-1);
            if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
        }
        defaultPath = util.getFullPathName(pathname);
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "setDefaultPath");
        exit(1);
    }
}
/*********************************************************************************************/
void CurrentFile::setTestFilePath(string pathname)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        if (pathname != "") {
            //add / to name if needed
            string lastChar = pathname.substr(pathname.length()-1);
            if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
        }
        
        testFilePath = util.getFullPathName(pathname);
        
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "setTestFilePath");
        exit(1);
    }
}
/*********************************************************************************************/
void CurrentFile::setBlastPath(string pathname)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        if (pathname != "") {
            //add / to name if needed
            string lastChar = pathname.substr(pathname.length()-1);
            if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
        }
        blastPath = util.getFullPathName(pathname);
        
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "setBlastPath");
        exit(1);
    }
}
/*********************************************************************************************/
void CurrentFile::setHomePath(string pathname)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        if (pathname != "") {
            //add / to name if needed
            string lastChar = pathname.substr(pathname.length()-1);
            if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
        }
        homePath = util.getFullPathName(pathname);
        m->setHomePath(homePath);
        
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "setHomePath");
        exit(1);
    }
}
/*********************************************************************************************/
void CurrentFile::setPaths(vector<string> pathVariables)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        for (int i = 0; i < pathVariables.size(); i++) {
            string pathname = pathVariables[i];
            if (pathname != "") {
                //add / to name if needed
                string lastChar = pathname.substr(pathname.length()-1);
                if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
            }
            pathVariables[i] = util.getFullPathName(pathname);
        }
        
        paths = pathVariables;
        m->setPaths(paths);
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "setPaths");
        exit(1);
    }
}
/*********************************************************************************************/
void CurrentFile::setToolsPath(string pathname)  {
    try {
        lock_guard<std::mutex> guard(currentProtector);
        
        if (pathname != "") {
            //add / to name if needed
            string lastChar = pathname.substr(pathname.length()-1);
            if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
        }
        toolsPath = util.getFullPathName(pathname);
        
    }
    catch(exception& e) {
        m->errorOut(e, "CurrentFile", "setToolsPath");
        exit(1);
    }
}


/*********************************************************************************************/
