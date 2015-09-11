###################################################
#
# Makefile for mothur
#
###################################################

#
# Macros
#

64BIT_VERSION ?= yes
USEREADLINE ?= yes
USECOMPRESSION ?= no
USEBOOST ?= yes
MOTHUR_FILES="\"Enter_your_default_path_here\""
RELEASE_DATE = "\"7/27/2015\""
VERSION = "\"1.36.1\""


ifeq  ($(strip $(64BIT_VERSION)),yes)
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
    LIBS += -lreadline
endif


#The boost libraries allow you to read gz files.
ifeq  ($(strip $(USEBOOST)),yes)
    LIBS += -lboost_iostreams
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

