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
MOTHUR_FILES="\"Enter_your_default_path_here\""
RELEASE_DATE = "\"3/31/2015\""
VERSION = "\"1.35.1\""
FORTAN_COMPILER = gfortran
FORTRAN_FLAGS =

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
    #CXXFLAGS += -mtune=native -march=native

    CXXFLAGS += -DBIT_VERSION
    FORTRAN_FLAGS = -m64
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

# if you want to enable reading and writing of compressed files, set to yes.
# The default is no.  this may only work on unix-like systems, not for windows.


ifeq  ($(strip $(USECOMPRESSION)),yes)
    CXXFLAGS += -DUSE_COMPRESSION
endif

#
# INCLUDE directories for mothur
#
#
    VPATH=calculators:chimera:classifier:clearcut:commands:communitytype:datastructures:metastats:randomforest:read:svm
    skipUchime := uchime_src/
    subdirs := $(filter-out  $(skipUchime), $(wildcard */))
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
    OBJECTS+=$(patsubst %.f,%.o,$(wildcard *.f))

mothur : fortranSource $(OBJECTS) uchime
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) -o $@ $(OBJECTS) $(LIBS)
	strip mothur

uchime:
	cd uchime_src && ./mk && mv uchime .. && cd ..

fortranSource:
	${FORTAN_COMPILER} -c $(FORTRAN_FLAGS) *.f

install : mothur uchime
	mv mothur ..
	mv uchime ..


%.o : %.c %.h
	$(COMPILE.c) $(OUTPUT_OPTION) $<
%.o : %.cpp %.h
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
%.o : %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<



clean :
	@rm -f $(OBJECTS)
	@rm -f uchime

