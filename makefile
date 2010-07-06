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

CYGWIN_BUILD ?= no
ifeq  ($(strip $(CYGWIN_BUILD)),yes)
    CXXFLAGS += -mno-cygwin
    LDFLAGS += -mno-cygwin 
#-lpsapi
endif

64BIT_VERSION ?= yes

ifeq  ($(strip $(64BIT_VERSION)),yes)
    TARGET_ARCH += -arch x86_64
endif

# if you do not want to use the readline library, set this to no.
# make sure you have the library installed

USEREADLINE ?= yes

ifeq  ($(strip $(USEREADLINE)),yes)
    CXXFLAGS += -DUSE_READLINE
    LDFLAGS += \
      -lreadline\
      -lncurses\
      -L../readline-6.0
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

clean :
	@rm -f $(OBJECTS)

