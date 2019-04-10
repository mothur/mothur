#include <stdlib.h>
#include <stdio.h>

#include <string.h>
#include "FileUtils.h"
#include "MatrixUtils.h"

#define SEPERATOR       " "
#define DELIM           ","
#define MAX_LINE_LENGTH 1024


/* simply removes trailing \r\n*/
char *strpret(char *szString)
{
  char *szTemp;

  szTemp = strpbrk(szString,"\n");

  if(szTemp!=NULL){
    *szTemp = '\0';
  }
  return szString;
}

char *  parseFromFile(const char *szFile, char *szParam)
{
  FILE *ifp    = NULL;
  char *szRet  = NULL;
  
  ifp = fopen(szFile,"r");
  
  if(ifp){
    char *szLine = NULL;
    char *szTok  = NULL;
    
    szLine = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
    if(szLine){
      
      while(fgets(szLine,MAX_LINE_LENGTH,ifp)){
	szLine = strpret(szLine);
	szTok = strtok(szLine,SEPERATOR);
	while(szTok){
	  if(strcmp(szTok,szParam)==0){
	    szTok = strtok(NULL,SEPERATOR);
	    if(strcmp(szTok,"=")==0){
	      szTok = strtok(NULL,SEPERATOR);
	      if(szTok){
		szRet = strdup(szTok);
		free(szLine);
		szLine = NULL;
		fclose(ifp);
		return szRet;
	      }
	    }
	
	  }
	
	  szTok = strtok(NULL,SEPERATOR);
	} /* end of while */
    
      } /* end of while */
    
      free(szLine);
      szLine = NULL;
    }
    else{
      fprintf(stderr,"Failed to allocate memory\n");
      exit(EXIT_FAILURE);
    }  
    
    fclose(ifp);
  }
  else{
    fprintf(stderr,"Failed to open file %s \n",szFile);
    exit(EXIT_FAILURE);
  }
  
  return szRet;
}

int extractFloatFromFile(const char *szFile,char *szName,float *pFloat)
{
  char *szTemp  = NULL;
  char *pcError = NULL;

  szTemp = parseFromFile(szFile,szName);
  if(szTemp == NULL)
    return 0;

  (*pFloat) = (float) strtod(szTemp,&pcError);
  if(*pcError!='\0'){
    free(szTemp);
    return 0;
  }
    
  free(szTemp);
  return 1;
}

int extractDoubleFromFile(const char *szFile,char *szName,double *pDouble)
{
  char *szTemp  = NULL;
  char *pcError = NULL;

  szTemp = parseFromFile(szFile,szName);
  if(szTemp == NULL)
    return 0;

  (*pDouble) = (double) strtod(szTemp,&pcError);
  if(*pcError!='\0'){
    free(szTemp);
    return 0;
  }
    
  free(szTemp);
  return 1;
}

int extractIntFromFile(const char *szFile,char *szName,int *pInt)
{
  char *szTemp  = NULL;
  char *pcError = NULL;

  szTemp = parseFromFile(szFile,szName);
  if(szTemp == NULL)
    return 0;
  
  (*pInt) = strtol(szTemp,&pcError,10);
  if(*pcError!='\0'){
    free(szTemp);
    return 0;
  }

  free(szTemp);
  return 1;
}

int extractLongFromFile(const char *szFile,char *szName,long *pLong)
{
  char *szTemp  = NULL;
  char *pcError = NULL;

  szTemp = parseFromFile(szFile,szName);
  if(szTemp == NULL)
    return 0;
  
  (*pLong) = strtol(szTemp,&pcError,10);
  if(*pcError!='\0'){
    free(szTemp);
    return 0;
  }

  free(szTemp);
  return 1;
}

int extractUShortsFromFile(const char *szFile, char *szName, unsigned short* aShorts, int nShorts)
{
  char *szTemp  = NULL;
  char *pcError = NULL;
  char *tok     = NULL;
  int  i        = 0;

  if(nShorts == 0)
    return 0;

  szTemp = parseFromFile(szFile,szName);
  if(szTemp == NULL)
    return 0;
  
  tok = strtok(szTemp, DELIM);
  
  aShorts[0] = (unsigned short) strtol(tok,&pcError,10);
  if(*pcError!='\0'){
    free(szTemp);
    return 0;
  }
  
  for(i = 1; i < nShorts; i++){
    tok = strtok(NULL,DELIM);
    aShorts[i] = (unsigned short) strtol(tok,&pcError,10);
    if(*pcError!='\0'){
      free(szTemp);
      return 0;
    }    
  }
  

  free(szTemp);
  return 1;
}

int extractDoublesFromFile(const char *szFile, char *szName, double* aDoubles, int nDoubles)
{
  char *szTemp  = NULL;
  char *pcError = NULL;
  char *tok     = NULL;
  int  i        = 0;

  if(nDoubles == 0)
    return 0;

  szTemp = parseFromFile(szFile,szName);
  if(szTemp == NULL)
    return 0;
  
  tok = strtok(szTemp, DELIM);
  
  aDoubles[0] =  strtod(tok,&pcError);
  if(*pcError!='\0'){
    free(szTemp);
    return 0;
  }
  
  for(i = 1; i < nDoubles; i++){
    tok = strtok(NULL,DELIM);
    aDoubles[i] = strtod(tok,&pcError);
    if(*pcError!='\0'){
      free(szTemp);
      return 0;
    }    
  }
  

  free(szTemp);
  return 1;
}

int extractIntMatrixFromFile(const char *szFile, char *szName, int** aanInts, int nRows, int nCols)
{
  FILE *ifp = fopen(szFile,"r");

  if(ifp){
    char szBuffer[MAX_LINE_LENGTH];
   
    while(fgets(szBuffer,MAX_LINE_LENGTH,ifp)){
      strpret(szBuffer);
      if(strcmp(szBuffer, szName) == 0){ /* i.e now found matrix*/
	if(!readMatrix2di(ifp,aanInts,nRows,nCols)){
	  fclose(ifp);
	  fprintf(stderr,"Failed to read 2d integer matrix %s in file %s\n",szName,szFile);
	  fprintf(stderr,"Lacked requisite number of rows and columns %d %d\n",nRows, nCols);
	  exit(EXIT_FAILURE);	  
	}
	else{
	  fclose(ifp);
	  return 1;
	}	  
      } /* name searching if*/
    } /*end line reading while*/
    
    fclose(ifp);
    fprintf(stderr,"Failed to find 2d integer matrix %s in file %s\n",szName,szFile);
    exit(EXIT_FAILURE);
  }
  else{
    fprintf(stderr,"Failed to open file %s for reading\n",szFile);
    exit(EXIT_FAILURE);
  }

}

int extractFloatMatrixFromFile(const char *szFile, char *szName, float** aafFloats, int nRows, int nCols)
{
  FILE *ifp = fopen(szFile,"r");

  if(ifp != NULL){
    char szBuffer[MAX_LINE_LENGTH];
    
    while(fgets(szBuffer,MAX_LINE_LENGTH,ifp)){
      strpret(szBuffer);
      if(strcmp(szBuffer, szName) == 0){ /* i.e now found matrix*/
	if(!readMatrix2df(ifp,aafFloats,nRows,nCols)){
	  fclose(ifp);
	  fprintf(stderr,"Failed to read 2d float matrix %s in file %s\n",szName,szFile);
	  fprintf(stderr,"Lacked requisite number of rows and columns %d %d\n",nRows, nCols);
	  exit(EXIT_FAILURE);	  
	}
	else{
	  fclose(ifp);
	  return 1;
	}	  
      } /* name searching if*/
    } /*end line reading while*/
    
    fclose(ifp);
    fprintf(stderr,"Failed to find 2d float matrix %s in file %s\n",szName,szFile);
    exit(EXIT_FAILURE);
  }
  else{
    fprintf(stderr,"Failed to open file %s for reading\n",szFile);
    exit(EXIT_FAILURE);
  }

}

int extractDoubleMatrixFromFile(const char *szFile, char *szName, double** aadDoubles, int nRows, int nCols)
{
  FILE *ifp = fopen(szFile,"r");

  if(ifp != NULL){
    char szBuffer[MAX_LINE_LENGTH];
    
    while(fgets(szBuffer,MAX_LINE_LENGTH,ifp)){
      strpret(szBuffer);
      if(strcmp(szBuffer, szName) == 0){ /* i.e now found matrix*/
	if(!readMatrix2dd(ifp,aadDoubles,nRows,nCols)){
	  fclose(ifp);
	  fprintf(stderr,"Failed to read 2d double matrix %s in file %s\n",szName,szFile);
	  fprintf(stderr,"Lacked requisite number of rows and columns %d %d\n",nRows, nCols);
	  exit(EXIT_FAILURE);	  
	}
	else{
	  fclose(ifp);
	  return 1;
	}	  
      } /* name searching if*/
    } /*end line reading while*/
    
    fclose(ifp);
    fprintf(stderr,"Failed to find 2d double matrix %s in file %s\n",szName,szFile);
    exit(EXIT_FAILURE);
  }
  else{
    fprintf(stderr,"Failed to open file %s for reading\n",szFile);
    exit(EXIT_FAILURE);
  }

}


