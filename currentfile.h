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
	
	
	void setListFile(string f)			{ listfile = f;				}
	void setTreeFile(string f)			{ treefile = f;				}
	void setGroupFile(string f)			{ groupfile = f;			}		
	void setPhylipFile(string f)		{ phylipfile = f;			}
	void setColumnFile(string f)		{ columnfile = f;			}
	void setNameFile(string f)			{ namefile = f;				}	
	void setRabundFile(string f)		{ rabundfile = f;			}
	void setSabundFile(string f)		{ sabundfile = f;			}
	void setSharedFile(string f)		{ sharedfile = f;			}
	void setRelAbundFile(string f)		{ relabundfile = f;			}
	void setOrderFile(string f)			{ orderfile = f;			}
	void setOrderGroupFile(string f)	{ ordergroupfile = f;		}
	void setDesignFile(string f)		{ designfile = f;			}
	void setFastaFile(string f)			{ fastafile = f;			}
	void setSFFFile(string f)			{ sfffile = f;				}
	void setQualFile(string f)			{ qualfile = f;				}
	void setOligosFile(string f)		{ oligosfile = f;			}
	
private:
	static CurrentFile* _uniqueInstance;
	CurrentFile( const CurrentFile& ); // Disable copy constructor
	void operator=( const CurrentFile& ); // Disable assignment operator
	CurrentFile();
	~CurrentFile();
	
	string phylipfile, columnfile, listfile, rabundfile, sabundfile, namefile, groupfile, designfile;
	string orderfile, treefile, sharedfile, ordergroupfile, relabundfile, fastafile, qualfile, sfffile, oligosfile;
	
	
};
/***********************************************/

#endif

