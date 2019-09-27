###################################################
#
# Makefile for mothur
#
###################################################

#
# Macros
#
# OPTIMIZE - yes will increase speed of executable.
# USEREADLINE - link with readline libraries.  Must have readline installed. Windows set to no.
# USEBOOST - link with boost libraries. Must install boost. Allows the make.contigs command to read .gz files.
# USEHDF5 - link with HDF5cpp libraries. Must install HDF5. Allows the biom.info command to read Biom format 2.0.
# USEGSL - link with GNU Scientific libraries. Must install GSL. Allows the estimiator.single command to find diversity estimates.
# HDF5_LIBRARY_DIR - location of HDF5 libraries
# HDF5_INCLUDE_DIR - location of HDF5 include files
# BOOST_LIBRARY_DIR - location of boost libraries
# BOOST_INCLUDE_DIR - location of boost include files
# GSL_LIBRARY_DIR - location of GSL libraries
# GSL_INCLUDE_DIR - location of GSL include files
# MOTHUR_FILES - The MOTHUR_FILES parameter is optional, but allows you to set a default location for mothur to look for input files it can't find. This is often used for reference files you want to store in one location separate from your data.

PREFIX := ${CURDIR} 

OPTIMIZE ?= yes
USEREADLINE ?= yes
USEBOOST ?= no
USEHDF5 ?= no
USEGSL ?= no
LOGFILE_NAME ?= no
BOOST_LIBRARY_DIR ?= "\"Enter_your_boost_library_path_here\""
BOOST_INCLUDE_DIR ?= "\"Enter_your_boost_include_path_here\""
HDF5_LIBRARY_DIR ?= "\"Enter_your_HDF5_library_path_here\""
HDF5_INCLUDE_DIR ?= "\"Enter_your_HDF5_include_path_here\""
GSL_LIBRARY_DIR ?= "\"Enter_your_GSL_library_path_here\""
GSL_INCLUDE_DIR ?= "\"Enter_your_GSL_include_path_here\""
MOTHUR_FILES="\"Enter_your_default_path_here\""
VERSION = "\"1.43.0\""


# Set a static logfile name
ifeq  ($(strip $(LOGFILE_NAME)),yes)
    LOGFILE_NAME="\"silent\""
endif

ifeq  ($(strip $(OPTIMIZE)),yes)
    CXXFLAGS += -O3
endif

CXXFLAGS += -std=c++11 -pthread -DVERSION=${VERSION}
LDFLAGS += -std=c++11 -pthread

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


#User specified boost library
ifeq  ($(strip $(USEBOOST)),yes)

    LDFLAGS += -L ${BOOST_LIBRARY_DIR}

    LIBS += -lboost_iostreams -lz
    CXXFLAGS += -DUSE_BOOST -I ${BOOST_INCLUDE_DIR}
endif

#User specified HDF5 library
ifeq  ($(strip $(USEHDF5)),yes)

LDFLAGS += -L ${HDF5_LIBRARY_DIR} -lhdf5 -lhdf5_cpp
CXXFLAGS += -DUSE_HDF5 -I ${HDF5_INCLUDE_DIR}

endif

#User specified GSL library
ifeq  ($(strip $(USEGSL)),yes)

LDFLAGS += -L ${GSL_LIBRARY_DIR} -lgsl -lgslcblas -lm
CXXFLAGS += -DUSE_GSL -I ${GSL_INCLUDE_DIR}

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
	cd source/uchime_src && export CXX=$(CXX) && ./mk && mv uchime ../../ && cd ..

install : mothur uchime
#if [ "${CURDIR}" = "$(PREFIX)" ]; then \
#		echo 'done'; \
#	else \
#		mkdir -p $(PREFIX); \
#		for file in mothur uchime; do \
#			cp -f $$file $(PREFIX) ; \
#		done \
#	fi


%.o : %.c %.h
	$(COMPILE.c) $(OUTPUT_OPTION) $<
%.o : %.cpp %.h
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
%.o : %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<



clean :
	@rm -f $(OBJECTS)
	@rm -f uchime
