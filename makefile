###################################################
#
# Makefile for mothur
# Created: June 29, 2010
#
###################################################

#
# Macros
#

# Optimize to level 3:

CXXFLAGS += -O3

64BIT_VERSION ?= yes

ifeq  ($(strip $(64BIT_VERSION)),yes)
	#if you are using centos uncomment the following lines
	#CXX = g++44
	
	#if you are a mac user use the following line
	TARGET_ARCH += -arch x86_64
	
	#if you are a linux user use the following line
	#CXXFLAGS += -mtune=native -march=native -m64
	
	 CXXFLAGS += -DBIT_VERSION
endif

MOTHUR_FILES = "\"../release\""

RELEASE_DATE = "\"10/26/2010\""
VERSION = "\"1.14.0\""

CXXFLAGS += -DRELEASE_DATE=${RELEASE_DATE} -DVERSION=${VERSION}

ifeq  ($(strip $(MOTHUR_FILES)),"\"Enter_your_default_path_here\"")
else
	CXXFLAGS += -DMOTHUR_FILES=${MOTHUR_FILES}
endif

CYGWIN_BUILD ?= no
ifeq  ($(strip $(CYGWIN_BUILD)),yes)
    CXXFLAGS += -mno-cygwin
    LDFLAGS += -mno-cygwin 
endif

# if you do not want to use the readline library, set this to no.
# make sure you have the library installed

USEREADLINE ?= yes

ifeq  ($(strip $(USEREADLINE)),yes)
    CXXFLAGS += -DUSE_READLINE
    LIBS = \
      -lreadline\
      -lncurses
endif

USEMPI ?= no

ifeq  ($(strip $(USEMPI)),yes)
    CXX = mpic++
    CXXFLAGS += -DUSE_MPI
endif

# if you want to enable reading and writing of compressed files, set to yes.
# The default is no.  this may only work on unix-like systems.

USECOMPRESSION ?= no

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

mothur : $(OBJECTS)
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) -o $@ $(OBJECTS) $(LIBS)

install : mothur
	cp mothur ../Release/mothur
	
%.o : %.c %.h
	$(COMPILE.c) $(OUTPUT_OPTION) $<
%.o : %.cpp %.h
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
%.o : %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<


clean :
	@rm -f $(OBJECTS)

