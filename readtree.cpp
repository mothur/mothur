/*
 *  readtree.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readtree.h"

/***********************************************************************/
//Parent Class
// The following functions are used by all reading formats.
/***********************************************************************/
ReadTree::ReadTree() { 
	try {
		globaldata = GlobalData::getInstance(); 
		T = new Tree(); 
		numNodes = T->getNumNodes();
		numLeaves = T->getNumLeaves();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function ReadTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function ReadTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************************/
int ReadTree::readSpecialChar(istream& f, char c, string name) {
    try {
		char d;
	
		while(isspace(d=f.get()))		{;}
		if(d == EOF){
			cerr << "Error: Input file ends prematurely, expecting a " << name << "\n";
			exit(1);
		}
		if(d != c){
			cerr << "Error: Expected " << name << " in input file.  Found " << d << ".\n";
			exit(1);
		}
		if(d == ')' && f.peek() == '\n'){
			while(isspace(d=f.get()))		{;}
			f.putback(d);
		}	
	
		return d;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function readSpecialChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function readSpecialChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/**************************************************************************************************/

int ReadTree::readNodeChar(istream& f) {
	try {
		char d;
		while(isspace(d=f.get()))		{;}
		if(d == EOF){
			cerr << "Error: Input file ends prematurely, expecting a left parenthesis\n";
			exit(1);
		}
		return d;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function readNodeChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function readNodeChar. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/**************************************************************************************************/

float ReadTree::readBranchLength(istream& f) {
    try {
		float b;
	
		if(!(f >> b)){
			cerr << "Error: Missing branch length in input tree.\n";
			exit(1);
		}
    
		return b;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function readBranchLength. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function readBranchLength. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}


/***********************************************************************/
/***********************************************************************/


//Child Classes Below

/***********************************************************************/
/***********************************************************************/
//This class reads a file in Newick form and stores it in a tree.

void ReadNewickTree::read() {
	try {
		int n = 0;
		int lc, rc; 
		
		int rooted = 0;
	
		int ch = filehandle.peek();	
		
		if(ch == '('){
			n = numLeaves;  //number of leaves / sequences, we want node 1 to start where the leaves left off
			lc = readNewickInt(filehandle, n, T);
		
			if(filehandle.peek()==','){							
				readSpecialChar(filehandle,',',"comma");
			}
			// ';' means end of tree.												
			else if((ch=filehandle.peek())==';' || ch=='['){		
				rooted = 1;									
			}												
			if(rooted != 1){								
				rc = readNewickInt(filehandle, n, T);
				if(filehandle.peek() == ')'){					
					readSpecialChar(filehandle,')',"right parenthesis");
				}											
			}												
		}
		//note: treeclimber had the code below added - not sure why?
		 else{
			filehandle.putback(ch);
			char name[MAX_LINE];
			filehandle.get(name, MAX_LINE,'\n');
			SKIPLINE(filehandle, ch);
		
		
			n = T->getIndex(name);
			if(n!=0){
				cerr << "Internal error: The only taxon is not taxon 0.\n";
				exit(1);
			}
			lc = rc = -1;
		} 
		
		while((ch=filehandle.get())!=';'){;}						
			if(rooted != 1){									
			T->tree[n].setChildren(lc,rc);
			T->tree[n].setBranchLength(0);
			T->tree[n].setParent(-1);
			if(lc!=-1){		T->tree[lc].setParent(n);		}
			if(rc!=-1){		T->tree[rc].setParent(n);		}
		}
		
		//save tree for later commands
		globaldata->gTree = T;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadNewickTree class Function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadNewickTree class function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/**************************************************************************************************/

int ReadNewickTree::readNewickInt(istream& f, int& n, Tree* T) {
	try {
		int c = readNodeChar(f);
    
		if(c == '('){
			int lc = readNewickInt(f, n, T);
			readSpecialChar(f,',',"comma");
		
			int rc = readNewickInt(f, n, T);		
			if(f.peek()==')'){	
				readSpecialChar(f,')',"right parenthesis");					
			}			
		
			if(f.peek() == ':'){									      
				readSpecialChar(f,':',"colon");							
				if(n >= numNodes){	cerr << "Error: Too many nodes in input tree\n";  exit(1); }
				T->tree[n].setBranchLength(readBranchLength(f));
			}else{T->tree[n].setBranchLength(0.0); }						
		
			T->tree[n].setChildren(lc,rc);
			T->tree[lc].setParent(n);
			T->tree[rc].setParent(n);
		
			return n++;
		}else{
			f.putback(c);
			string name = "";
			char d=f.get();
			while(d != ':' && d != ',' && d!=')' && d!='\n'){					
				name += d;
				d=f.get();
			}
		
			int blen = 0;
			if(d == ':')	{		blen = 1;			}		
		
			f.putback(d);
		
			//set group info
			string group = globaldata->gTreemap->getGroup(name);
			
			//find index in tree of name
			int n1 = T->getIndex(name);
			
			if(n1 == -1){cerr << "Name: " << name << " not found\n"; exit(1);}
			
			else T->tree[n1].setGroup(group);
		
			T->tree[n1].setChildren(-1,-1);
		
			if(blen == 1){	
				f.get();		
				T->tree[n1].setBranchLength(readBranchLength(f));
			}else{
				T->tree[n1].setBranchLength(0.0);
			}
		
			while((c=f.get())!=0 && (c != ':' && c != ',' && c!=')') )		{;}		
			f.putback(c);
		
			return n1;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadNewickTree class Function readNewickInt. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadNewickTree class function readNewickInt. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/**************************************************************************************************/
/**************************************************************************************************/

