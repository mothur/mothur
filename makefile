###################################################
#
# Makefile for mothur
# Created: June 29, 2010
#
###################################################

#
# Macros
#

USEMPI ?= no
64BIT_VERSION ?= yes
USEREADLINE ?= no
CYGWIN_BUILD ?= no
USECOMPRESSION ?= no
MOTHUR_FILES="\"Enter_your_default_path_here\""
RELEASE_DATE = "\"3/13/2012\""
VERSION = "\"1.24.0\""
FORTAN_COMPILER = gfortran

# Optimize to level 3:
CXXFLAGS += -O3 

ifeq  ($(strip $(64BIT_VERSION)),yes)
	#if you are using centos uncomment the following lines
	#CXX = g++44
	
	#if you are a mac user use the following line
	TARGET_ARCH += -arch x86_64
	
	#if you using cygwin to build Windows the following line
	#CXX = x86_64-w64-mingw32-g++
	#CC = x86_64-w64-mingw32-g++
	#FORTAN_COMPILER = x86_64-w64-mingw32-gfortran
 	#TARGET_ARCH += -m64 -static

	#if you are a linux user use the following line
	#CXXFLAGS += -mtune=native -march=native -m64
	
	CXXFLAGS += -DBIT_VERSION
endif


CXXFLAGS += -DRELEASE_DATE=${RELEASE_DATE} -DVERSION=${VERSION}

ifeq  ($(strip $(MOTHUR_FILES)),"\"Enter_your_default_path_here\"")
else
	CXXFLAGS += -DMOTHUR_FILES=${MOTHUR_FILES}
endif

ifeq  ($(strip $(CYGWIN_BUILD)),yes)
    CXXFLAGS += -mno-cygwin
    LDFLAGS += -mno-cygwin 
endif

# if you do not want to use the readline library, set this to no.
# make sure you have the library installed


ifeq  ($(strip $(USEREADLINE)),yes)
    CXXFLAGS += -DUSE_READLINE
    LIBS = \
      -lreadline\
      -lncurses
endif


ifeq  ($(strip $(USEMPI)),yes)
    CXX = mpic++
    CXXFLAGS += -DUSE_MPI
endif

# if you want to enable reading and writing of compressed files, set to yes.
# The default is no.  this may only work on unix-like systems, not for windows.


ifeq  ($(strip $(USECOMPRESSION)),yes)
  CXXFLAGS += -DUSE_COMPRESSION
endif

#
# INCLUDE directories for mothur
#

     CXXFLAGS += -I.

#
# Get the list of all .cpp files, rename to .o files
#

OBJECTS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
OBJECTS+=$(patsubst %.c,%.o,$(wildcard *.c))
OBJECTS+=$(patsubst %.f,%.o,$(wildcard *.f))

mothur : fortranSource $(OBJECTS) uchime
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) -o $@ $(OBJECTS) $(LIBS)
	
	strip mothur
 
uchime:
	cd uchime_src && ./mk && mv uchime .. && cd ..
	
fortranSource:
	${FORTAN_COMPILER} -c -m64 *.f

install : mothur
#	cp mothur ../Release/mothur
	
%.o : %.c %.h
	$(COMPILE.c) $(OUTPUT_OPTION) $<
%.o : %.cpp %.h
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
%.o : %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<


clean :
	@rm -f $(OBJECTS)
	@rm -f uchime

