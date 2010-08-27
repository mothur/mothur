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

MOTHUR_FILES = "\"../Release\""

RELEASE_DATE = "\"8/30/2010\""
VERSION = "\"1.13.0\""

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

64BIT_VERSION ?= yes

ifeq  ($(strip $(64BIT_VERSION)),yes)
    TARGET_ARCH += -arch x86_64
	 CXXFLAGS += -DBIT_VERSION
	
	#if you are using centos uncomment the following lines
	#CXX = g++44
	#CXXFLAGS += -mtune=native -march=native -m64
endif

# if you do not want to use the readline library, set this to no.
# make sure you have the library installed

USEREADLINE ?= yes

ifeq  ($(strip $(USEREADLINE)),yes)
    CXXFLAGS += -DUSE_READLINE
    LDFLAGS += \
      -lreadline\
      -lncurses
endif

USEMPI ?= no

ifeq  ($(strip $(USEMPI)),yes)
    CXX = mpic++
    CXXFLAGS += -DUSE_MPI
endif

#
# INCLUDE directories for mothur
#

     CXXFLAGS += -I.

#
# Get the list of all .cpp files, rename to .o files
#

OBJECTS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

mothur : $(OBJECTS)
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) -o $@ $(OBJECTS)

install : mothur
	cp mothur ../Release/mothur

%.o : %.cpp %.h
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
%.o : %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<


clean :
	@rm -f $(OBJECTS)

