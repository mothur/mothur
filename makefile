###################################################
#
# Makefile for mothur
#
###################################################

#
# Macros
#

USEMPI ?= no
64BIT_VERSION ?= yes
USEREADLINE ?= yes
USECOMPRESSION ?= no
USEBOOST ?= yes
MOTHUR_FILES="\"Enter_your_default_path_here\""
RELEASE_DATE = "\"7/23/2015\""
VERSION = "\"1.36.0\""

# Optimize to level 3:
CXXFLAGS += -O3


ifeq  ($(strip $(64BIT_VERSION)),yes)
    #if you are a mac user use the following line
    TARGET_ARCH += -arch x86_64

    #if you using cygwin to build Windows the following line
    #CXX = x86_64-w64-mingw32-g++
    #CC = x86_64-w64-mingw32-g++
    #TARGET_ARCH += -m64 -static

    #if you are a linux user use the following line
    #CXXFLAGS += -mtune=native -march=native

    CXXFLAGS += -DBIT_VERSION
endif


CXXFLAGS += -DRELEASE_DATE=${RELEASE_DATE} -DVERSION=${VERSION}

ifeq  ($(strip $(MOTHUR_FILES)),"\"Enter_your_default_path_here\"")
else
    CXXFLAGS += -DMOTHUR_FILES=${MOTHUR_FILES}
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

#The boost libraries allow you to read gz files.
ifeq  ($(strip $(USEBOOST)),yes)
    BOOST_INCLUDE_DIR="/usr/local/include"
    BOOST_LIBRARY_DIR="/usr/local/lib"

    CXXFLAGS += -DUSE_BOOST

    LIBS += \
    ${BOOST_LIBRARY_DIR}/libboost_iostreams.a \
    ${BOOST_LIBRARY_DIR}/zlib.a

    #if linux or windows then ${BOOST_LIBRARY_DIR}/libz.a
endif


# if you want to enable reading and writing of compressed files, set to yes.
# The default is no.  this may only work on unix-like systems, not for windows.


ifeq  ($(strip $(USECOMPRESSION)),yes)
    CXXFLAGS += -DUSE_COMPRESSION
endif

#
# INCLUDE directories for mothur
#
#
    VPATH=source/calculators:source/chimera:source/classifier:source/clearcut:source/commands:source/communitytype:source/datastructures:source/metastats:source/randomforest:source/read:source/svm
    skipUchime := source/uchime_src/
    subdirs :=  $(sort $(dir $(filter-out  $(skipUchime), $(wildcard source/*/))))
    subDirIncludes = $(patsubst %, -I %, $(subdirs))
    subDirLinking =  $(patsubst %, -L%, $(subdirs))
    CXXFLAGS += -I. $(subDirIncludes)
    LDFLAGS += $(subDirLinking)


#
# Get the list of all .cpp files, rename to .o files
#
    OBJECTS=$(patsubst %.cpp,%.o,$(wildcard $(addsuffix *.cpp,$(subdirs))))
    OBJECTS+=$(patsubst %.c,%.o,$(wildcard $(addsuffix *.c,$(subdirs))))
    OBJECTS+=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
    OBJECTS+=$(patsubst %.c,%.o,$(wildcard *.c))

mothur : $(OBJECTS) uchime
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) -o $@ $(OBJECTS) $(LIBS)
	strip mothur


uchime:
	cd source/uchime_src && ./mk && mv uchime ../../ && cd ..


install : mothur


%.o : %.c %.h
	$(COMPILE.c) $(OUTPUT_OPTION) $<
%.o : %.cpp %.h
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
%.o : %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<



clean :
	@rm -f $(OBJECTS)
	@rm -f uchime

