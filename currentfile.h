#ifndef CURRENTFILE_H
#define CURRENTFILE_H

/*
 *  currentfile.h
 *  Mothur
 *
 *  Created by westcott on 3/15/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"
#include "mothur.h"

/***********************************************/

class CurrentFile {
	
public:
	
	static CurrentFile* getInstance();
	
	string getPhylipFile()		{ return phylipfile;		}
	string getColumnFile()		{ return columnfile;		}
	string getListFile()		{ return listfile;			}
	string getRabundFile()		{ return rabundfile;		}
	string getSabundFile()		{ return sabundfile;		}
	string getNameFile()		{ return namefile;			}	
	string getGroupFile()		{ return groupfile;			}	
	string getOrderFile()		{ return orderfile;			}
	string getOrderGroupFile()	{ return ordergroupfile;	}
	string getTreeFile()		{ return treefile;			}
	string getSharedFile()		{ return sharedfile;		}
	string getRelAbundFile()	{ return relabundfile;		}
	string getDesignFile()		{ return designfile;		}
	string getFastaFile()		{ return fastafile;			}
	string getSFFFile()			{ return sfffile;			}
	string getQualFile()		{ return qualfile;			}
	string getOligosFile()		{ return oligosfile;		}
	
	
	void setListFile(string f)			{ listfile = m->getFullPathName(f);				}
	void setTreeFile(string f)			{ treefile = m->getFullPathName(f);				}
	void setGroupFile(string f)			{ groupfile = m->getFullPathName(f);			}		
	void setPhylipFile(string f)		{ phylipfile = m->getFullPathName(f);			}
	void setColumnFile(string f)		{ columnfile = m->getFullPathName(f);			}
	void setNameFile(string f)			{ namefile = m->getFullPathName(f);				}	
	void setRabundFile(string f)		{ rabundfile = m->getFullPathName(f);			}
	void setSabundFile(string f)		{ sabundfile = m->getFullPathName(f);			}
	void setSharedFile(string f)		{ sharedfile = m->getFullPathName(f);			}
	void setRelAbundFile(string f)		{ relabundfile = m->getFullPathName(f);			}
	void setOrderFile(string f)			{ orderfile = m->getFullPathName(f);			}
	void setOrderGroupFile(string f)	{ ordergroupfile = m->getFullPathName(f);		}
	void setDesignFile(string f)		{ designfile = m->getFullPathName(f);			}
	void setFastaFile(string f)			{ fastafile = m->getFullPathName(f);			}
	void setSFFFile(string f)			{ sfffile = m->getFullPathName(f);				}
	void setQualFile(string f)			{ qualfile = m->getFullPathName(f);				}
	void setOligosFile(string f)		{ oligosfile = m->getFullPathName(f);			}
	
private:
	
	static CurrentFile* _uniqueInstance;
	CurrentFile( const CurrentFile& ); // Disable copy constructor
	void operator=( const CurrentFile& ); // Disable assignment operator
	CurrentFile();
	~CurrentFile();
	
	MothurOut* m;
	string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile;
	string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile;
	
	
};
/***********************************************/

#endif

